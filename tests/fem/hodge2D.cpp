

/**/

#include <cmath>

#include <functional>
#include <vector>
#include <string>

#include "../../base/include.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.veewedgehodge.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: (2D) Hodge star, wedge, vee" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    auto M = UnitSquare2D();
    
    M.check();
    
    const int n = M.getinnerdimension();
    
    


    std::vector<std::vector<std::function<FloatVector(const FloatVector&)>>> experiments_field(n+1);
    std::vector<std::vector<Float>>                                          experiments_value(n+1);
    
    experiments_field[0].push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            return FloatVector({ 
                1.
            });
        }
    );
    experiments_value[0].push_back( 1. );

    
    experiments_field[0].push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            return FloatVector({ 
                vec[0]*vec[1]*vec[1] + std::exp( vec[0] * vec[1] )
            });
        }
    );
    experiments_value[0].push_back( 1./3. * 1./5. + square( Constants::euler*Constants::euler/2. - 1. ));

    
    experiments_field[1].push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            return FloatVector({ 
                1., 1.
            });
        }
    );
    experiments_value[1].push_back( sqrt(2.) );

    experiments_field[1].push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            return FloatVector({ 
                vec[0]*vec[1]*vec[1] + std::exp( vec[0] * vec[1] ),
                5.
            });
        }
    );
    experiments_value[1].push_back( 1./3. * 1./5. + square( Constants::euler*Constants::euler/2. - 1. ) + 25. );


    
    experiments_field[2] = experiments_field[0];
    experiments_value[2] = experiments_value[0];
    
    
    
    
    
    
    
    
    
    const int r_min = 0;
    
    const int r_max = 2;
    
    const int l_min = 0;
    
    const int l_max = 2;

    const int number_of_samples = 1;
    
        
    
    Float errors_mass [ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float errors_vee  [ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float errors_wedge[ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    
    
        
    for( int l = 0; l < l_min; l++ )
        M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ ){
        
        LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
        
        for( int k = 0; k <= n; k++ ) 
        for( int r = r_min; r <= r_max; r++ ) 
        {
            
            LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;

            LOG << "Form degree: " << space << k << nl;

            LOG << "assemble matrices..." << nl;
    
            errors_mass [k][ l-l_min ][ r-r_min ] = 0.;
            errors_vee  [k][ l-l_min ][ r-r_min ] = 0.;
            errors_wedge[k][ l-l_min ][ r-r_min ] = 0.;
            
            SparseMatrix broken_mass_primal_matrix = FEECBrokenMassMatrix ( M, M.getinnerdimension(),   k, r );
            SparseMatrix broken_mass_dual_matrix   = FEECBrokenMassMatrix ( M, M.getinnerdimension(), n-k, r );
            SparseMatrix broken_mass_scalar_matrix = FEECBrokenMassMatrix ( M, M.getinnerdimension(),   0, 2*r );
            SparseMatrix broken_mass_volume_matrix = FEECBrokenMassMatrix ( M, M.getinnerdimension(),   n, 2*r );

            FloatVector volume_integrals = FEECVolumeFormIntegral( M, M.getinnerdimension(), 2*r );

            FloatVector scalar_integrals = FEECScalarIntegral( M, M.getinnerdimension(), 2*r );

            SparseMatrix broken_hodge_primal_matrix  = FEECBrokenHodgeMatrix( M, M.getinnerdimension(), k, r );
            
            assert( broken_hodge_primal_matrix.getdimin()  == broken_mass_primal_matrix.getdimin() );
            assert( broken_hodge_primal_matrix.getdimout() == broken_mass_dual_matrix.getdimout()  );

            LOG << "sampling..." << nl;
    
            for( int i = 0; i < experiments_field[k].size(); i++ ){

                // LOG << "compare mass of u and star u...\t";
                
                // auto field = broken_mass_primal_matrix.createinputvector();
                // field.random();
                // field.normalize();

                auto current_field = experiments_field[k][i];
                auto current_value = experiments_value[k][i];

                auto field = Interpolation( M, n, k, r, current_field );
                
                assert( field.is_finite() );

                const auto dual_field = broken_hodge_primal_matrix * field;
                
                assert( dual_field.is_finite() );

                const Float mass_primal = field * ( broken_mass_primal_matrix * field );
                
                const Float mass_dual   = dual_field * ( broken_mass_dual_matrix * dual_field   );
                
                const auto error_mass = absolute( mass_primal - mass_dual );

                assert( std::isfinite(error_mass) );

                // LOG << mass_primal << space << mass_dual << space << mass_dual / mass_primal << space << error_mass << nl;
                
                errors_mass[k][l-l_min][r-r_min] = maximum( errors_mass[k][l-l_min][r-r_min], error_mass );

                {
                    // LOG << "integrate u vee u\t\t";

                    SparseMatrix broken_vee_with_primal_matrix  = FEECBrokenVeeMatrix( M, M.getinnerdimension(), k, r, k, r, field );

                    const auto vee_scalar = broken_vee_with_primal_matrix * field;

                    const Float vee_mass = scalar_integrals * vee_scalar;
                    // const Float vee_mass = vee_scalar * ( broken_mass_scalar_matrix * vee_scalar );

                    // LOG << mass_primal << space << vee_mass << space << vee_mass / mass_primal << nl;

                    Float error_vee = absolute( vee_mass - mass_primal );
                    
                    errors_vee[k][l-l_min][r-r_min] = maximum( errors_vee[k][l-l_min][r-r_min], error_vee );
                }

                {
                    LOG << "integrate u wedge star u\t"; 

                    SparseMatrix broken_wedge_with_primal_matrix  = FEECBrokenWedgeMatrix( M, M.getinnerdimension(), n-k, r, k, r, field );

                    const auto wedge_volume = broken_wedge_with_primal_matrix * dual_field;

                    Assert( volume_integrals.getdimension() == wedge_volume.getdimension(), volume_integrals.getdimension(), wedge_volume.getdimension() );
                    const Float wedge_mass = volume_integrals * wedge_volume;
                    // const Float wedge_mass = wedge_volume * ( broken_mass_volume_matrix * wedge_volume );

                    // if( k == n )
                    LOGPRINTF("primal mass %le vs wedge mass %le, ratio=%le\n", mass_primal, wedge_mass, wedge_mass / mass_primal );

                    Float error_wedge = absolute( wedge_mass - mass_primal );
                    
                    errors_wedge[k][l-l_min][r-r_min] = maximum( errors_wedge[k][l-l_min][r-r_min], error_wedge );
                }

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

    std::vector<ConvergenceTable> contables_mass(n + 1);
    std::vector<ConvergenceTable> contables_vee(n + 1);
    std::vector<ConvergenceTable> contables_wedge(n + 1);

    for (int k = 0; k <= n; k++) {
        contables_mass[k].table_name  = "Rounding errors D2K_mass" + std::to_string(k);
        contables_vee[k].table_name   = "Rounding errors D2K_vee" + std::to_string(k);
        contables_wedge[k].table_name = "Rounding errors D2K_wedge" + std::to_string(k);
    }

    for (int k = 0; k <= n; k++) {
        for (int r = r_min; r <= r_max; r++) {
            contables_mass[k]  << ("R" + std::to_string(r));
            contables_vee[k]   << ("R" + std::to_string(r));
            contables_wedge[k] << ("R" + std::to_string(r));
        }
        contables_mass[k]  << nl;
        contables_vee[k]   << nl;
        contables_wedge[k] << nl;
    }

    for (int k = 0; k <= n; k++)
    for (int l = l_min; l <= l_max; l++) 
    {
        for (int r = r_min; r <= r_max; r++) {
            contables_mass[k]  << errors_mass [k][l - l_min][r - r_min] / machine_epsilon;
            contables_vee[k]   << errors_vee  [k][l - l_min][r - r_min] / machine_epsilon;
            contables_wedge[k] << errors_wedge[k][l - l_min][r - r_min] / machine_epsilon;
        }
        contables_mass[k]  << nl;
        contables_vee[k]   << nl;
        contables_wedge[k] << nl;
    }

    for (int k = 0; k <= n; k++) {
        contables_mass[k].lg();
        contables_vee[k].lg();
        contables_wedge[k].lg();
        LOG << "                   " << nl;
    }

    
    
    LOG << "Check that differences are below: " << desired_closeness << nl;
    
    for( int l = l_min; l <= l_max; l++ ) 
    for( int r = r_min; r <= r_max; r++ ) 
    for( int k = 0; k <= n; k++ ) 
    {
        Assert( errors_mass[k][l-l_min][r-r_min] < desired_closeness, errors_mass[k][l-l_min][r-r_min], desired_closeness );
        Assert( errors_vee[k][l-l_min][r-r_min] < desired_closeness, errors_vee[k][l-l_min][r-r_min], desired_closeness );
        Assert( errors_wedge[k][l-l_min][r-r_min] < desired_closeness, errors_wedge[k][l-l_min][r-r_min], desired_closeness );
    }
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
