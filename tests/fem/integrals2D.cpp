

/**/

#include <cmath>

#include <functional>
#include <vector>
#include <string>

#include "../../base/include.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../fem/global.veewedgehodge.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: (2D) Auxiliary functions for integrating scalar fields and volume forms" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    // auto M = UnitSquare2D_simple();
    auto M = UnitSquare2D_strange14();
    // auto M = UnitTriangle2D();
    
    M.check();

    assert( M.getinnerdimension() == M.getouterdimension() );
    




    std::function<FloatVector(const FloatVector&)> constant_one
            = [](const FloatVector& vec) -> FloatVector {
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 1. });
                };
        
    
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
    std::vector<Float>                                          experiments_scalar_value;
    
    experiments_scalar_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            return FloatVector({ 
                1.
            });
        }
    );
    experiments_scalar_value.push_back( 1. );

    experiments_scalar_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            return FloatVector({ 
                std::sin( Constants::twopi * vec[0] ) * std::sin( Constants::twopi * vec[1] ) 
            });
        }
    );
    experiments_scalar_value.push_back( 0. );

    experiments_scalar_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            return FloatVector({ std::exp( vec[0] + vec[1] ) });
        }
    );
    experiments_scalar_value.push_back( std::exp(2.) - 2 * std::exp(1.) + 1. );
    experiments_scalar_value.push_back( 1. );

    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field = experiments_scalar_field;
    std::vector<Float>                                          experiments_volume_value = experiments_scalar_value;




    const int r_min = 0;
    
    const int r_max = 3;
    
    const int l_min = 0;
    
    const int l_max = 5;

    // const int n = 2;
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

            errors_scalar[ l-l_min ][ r-r_min ] = 0.;
            errors_volume[ l-l_min ][ r-r_min ] = 0.;
            
            FloatVector scalar_integrals = FEECScalarIntegral( M, M.getinnerdimension(), r );

            FloatVector volume_integrals = FEECVolumeFormIntegral( M, M.getinnerdimension(), r );

            // 1. check that the integral of scalar 1 is correct
            {
                FloatVector interpol_scalar_one  = Interpolation( M, M.getinnerdimension(), 0, r, constant_one );
            
                Float value = scalar_integrals * interpol_scalar_one;

                Assert( is_numerically_close( value, 1. ), value );
            }

            // 2. check that the integral of volume 1 is correct
            {
                FloatVector interpol_volume_one  = Interpolation( M, M.getinnerdimension(), n, r, constant_one );
            
                Float value = volume_integrals * interpol_volume_one;
                
                Assert( is_numerically_close( value, 1. ), value );
            }

            // 3. Find the integrals of the scalar fields 
            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
            {    
                auto scalarfield = experiments_scalar_field[i];

                auto interpol = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                auto interpol_integral = scalar_integrals * interpol;

                auto error = absolute( interpol_integral - experiments_scalar_value[i] );
                
                assert( std::isfinite(error) );

                errors_scalar[l-l_min][r-r_min] = maximum( errors_scalar[l-l_min][r-r_min], error );   
            }

            // 4. Find the integrals of the volume fields 
            
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
            {    
                auto volumefield = experiments_volume_field[i];

                auto interpol = Interpolation( M, M.getinnerdimension(), n, r, volumefield );

                auto interpol_integral = volume_integrals * interpol;

                auto error = absolute( interpol_integral - experiments_volume_value[i] );
                
                assert( std::isfinite(error) );

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
        if( l < l_max or r < r_max ) continue;

        Assert( errors_scalar[l-l_min][r-r_min] < desired_closeness, errors_scalar[l-l_min][r-r_min], desired_closeness );
        Assert( errors_volume[l-l_min][r-r_min] < desired_closeness, errors_volume[l-l_min][r-r_min], desired_closeness );
    }
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
