

/**/

#include <cmath>

#include <functional>
#include <vector>
#include <string>

#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/global.veewedgehodge.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    
    LOG << "Unit Test: (3D) Integrating scalar fields and volume forms" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    auto M = UnitCube3D();
    
    M.check();

    assert( M.getinnerdimension() == M.getouterdimension() );
    




    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
    std::vector<Float>                                          experiments_scalar_value;
    
    experiments_scalar_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 3 );
            return FloatVector({ 
                1.
            });
        }
    );
    experiments_scalar_value.push_back( 1. );

    experiments_scalar_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 3 );
            return FloatVector({ 
                std::sin( Constants::twopi * vec[0] ) * std::sin( Constants::twopi * vec[1] ) * std::sin( Constants::twopi * vec[2] ) 
            });
        }
    );
    experiments_scalar_value.push_back( 0. );

    experiments_scalar_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 3 );
            return FloatVector({ std::exp( vec[0] + vec[1] + vec[2] ) });
        }
    );
    experiments_scalar_value.push_back( power_numerical( Constants::euler - 1., 3 ) );
    
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field = experiments_scalar_field;
    std::vector<Float>                                          experiments_volume_value = experiments_scalar_value;




    const int r_min = 0;
    
    const int r_max = 2;
    
    const int l_min = 0;
    
    const int l_max = 2;

    // const int n = 3;
    const int n = M.getinnerdimension();
    assert( n == M.getinnerdimension() );
    
    Float errors_volume[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float errors_scalar[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    
    
        
    for( int l = 0; l < l_min; l++ )
        M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ ){
        
        LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
        
        for( int r = r_min; r <= r_max; r++ ) 
        {
            
            LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;

            LOG << "assemble matrices..." << nl;
    
            errors_scalar[ l-l_min ][ r-r_min ] = 0.;
            errors_volume[ l-l_min ][ r-r_min ] = 0.;
            
            FloatVector scalar_integrals = FEECScalarIntegral( M, M.getinnerdimension(), r );

            FloatVector volume_integrals = FEECVolumeFormIntegral( M, M.getinnerdimension(), r );

            if( r <= 1 )
            {
                FloatVector unitvector = FloatVector( scalar_integrals.getdimension(), 1.);
                Float unit = scalar_integrals * unitvector;
                Assert( is_numerically_close( unit, 1. ), unit );
            }

            if( r <= 1 )
            {
                FloatVector unitvector = FloatVector( volume_integrals.getdimension(), 1. );
                
                Float value = volume_integrals * unitvector;
                
                int count_positive = 0;
                for( int s = 0; s < M.count_simplices(n); s++ ) {
                    if( M.getOrientation(s) == 1. ) 
                        count_positive++;
                }
                Float desired_value_1 = count_positive - ( M.count_simplices(n) - count_positive );

                Float desired_value_2 = ( 1 - (n % 2) ) / factorial_numerical(n);
                
                Float desired_value = desired_value_1 * desired_value_2;
                
                Assert( is_numerically_close( value, desired_value ), 
                        value, desired_value, 
                        volume_integrals,
                        M.count_simplices(n) );
            }

            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
            {    
                auto scalarfield = experiments_scalar_field[i];

                auto interpol = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                auto interpol_integral = scalar_integrals * interpol;

                auto error = absolute( interpol_integral - experiments_scalar_value[i] );
                
                assert( isfinite(error) );

                errors_scalar[l-l_min][r-r_min] = maximum( errors_scalar[l-l_min][r-r_min], error );   
            }

            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
            {    
                auto volumefield = experiments_volume_field[i];

                auto interpol = Interpolation( M, M.getinnerdimension(), n, r, volumefield );

                // LOG << n << space << volume_integrals.getdimension() << space << interpol.getdimension() << nl;
                auto interpol_integral = volume_integrals * interpol;

                // LOG << interpol_integral << space << experiments_volume_value[i] << nl;
                // if( r == 0 and l == 0 ) { LOG << volume_integrals << nl << interpol << nl; }

                auto error = absolute( interpol_integral - experiments_volume_value[i] );
                
                assert( isfinite(error) );

                errors_volume[l-l_min][r-r_min] = maximum( errors_volume[l-l_min][r-r_min], error );   
            }
            
        }

        if( l != l_max )
        {
            LOG << "Refinement..." << nl;
        
            M.uniformrefinement();

            M.shake_interior_vertices();
        }
        
    } 
    
    
    
    LOG << "Convergence tables" << nl;

    ConvergenceTable contables_scalar;
    ConvergenceTable contables_volume;

    contables_scalar.table_name = "Rounding errors scalar";
    contables_volume.table_name = "Rounding errors volume";
    
    for (int r = r_min; r <= r_max; r++) {
        contables_scalar << ("R" + std::to_string(r));
        contables_volume << ("R" + std::to_string(r));
    }
    contables_scalar << nl;
    contables_volume << nl;
    
    for (int l = l_min; l <= l_max; l++) 
    {
        for (int r = r_min; r <= r_max; r++) {
            contables_scalar << errors_scalar[l - l_min][r - r_min];
            contables_volume << errors_volume[l - l_min][r - r_min];
        }
        contables_scalar << nl;
        contables_volume << nl;
    }

    contables_scalar.lg();
    contables_volume.lg();
    LOG << "                   " << nl;
    
    LOG << "Check that differences are below: " << desired_closeness << nl;
    
    for( int l = l_min; l <= l_max; l++ ) 
    for( int r = r_min; r <= r_max; r++ ) 
    {
        Assert( errors_scalar[l-l_min][r-r_min] < desired_closeness, errors_scalar[l-l_min][r-r_min], desired_closeness );
        Assert( errors_volume[l-l_min][r-r_min] < desired_closeness, errors_volume[l-l_min][r-r_min], desired_closeness );
    }
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
