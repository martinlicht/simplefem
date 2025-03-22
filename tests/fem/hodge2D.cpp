

/**/

#include <cmath>

#include <functional>
#include <vector>
#include <string>

#include "../../base/include.hpp"
#include "../../operators/composedoperators.hpp"
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
    M.getCoordinates().scale(1.5);
    M.check();
    
    const int n = M.getinnerdimension();
    
    

    /* 0. Before anything else, we test the Euclidean Hodge matrix */
    
    for( int d = 1; d <= 5; d++ )
    for( int k = 0; k <= d; k++ )
    {
        auto hodge = EuclideanHodgeStar( d,   k );
        auto other = EuclideanHodgeStar( d, d-k );

        auto prod = other * hodge;

        auto signed_identity = sign_power( k * (d-k) ) * IdentityMatrix( other.getdimout() );

        LOG << prod << nl << signed_identity << nl;
        Assert( ( prod - signed_identity ).norm() == 0., prod, signed_identity );
    }
    
    
    
    
    /* Prepare some scalar/vector/pseudo/volume fields with known mass */

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
    
    const int r_max = 3;
    
    const int l_min = 0;
    
    const int l_max = 3;

    const int number_of_random_trial_samples = 5;
    
    const int number_of_random_test_samples  = 5;
    
    Float errors_mass_dual [ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float errors_vee       [ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float errors_wedgeA    [ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float errors_wedgeB    [ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    
    Float sampled_mass_errors[ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    
    
    
        
    for( int l = 0; l < l_min; l++ )
        M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ ){
        
        for( int k = 0; k <= n; k++ ) 
        for( int r = r_min; r <= r_max; r++ ) 
        {
            
            LOG << "Level: " << l_min << " <= " << l << " <= " << l_max << nl;
        
            LOG << "Polydegree: " << r_min << " <= " << r << " <= " << r_max << nl;

            LOG << "Form degree: " << k << nl;

            LOG << "assemble matrices..." << nl;
            
            SparseMatrix broken_mass_primal_op = FEECBrokenMassMatrix ( M, M.getinnerdimension(),   k, r );
            SparseMatrix broken_mass_dual_op   = FEECBrokenMassMatrix ( M, M.getinnerdimension(), n-k, r );
            
            SparseMatrix broken_hodge_primal_op  = FEECBrokenHodgeStarMatrix( M, M.getinnerdimension(),   k, r );
            SparseMatrix broken_hodge_dual_op    = FEECBrokenHodgeStarMatrix( M, M.getinnerdimension(), n-k, r );
            
            SparseMatrix& alternative_hodge_star  = broken_hodge_primal_op; 
            // SparseMatrix alternative_hodge_star  = FEECBrokenHodgeStarMatrix( M, M.getinnerdimension(), k, r ); // TODO: Fix
            
            FloatVector volume_integrals = FEECVolumeFormIntegral( M, M.getinnerdimension(), 2*r );
            FloatVector scalar_integrals = FEECScalarIntegral(     M, M.getinnerdimension(), 2*r );

            assert( broken_hodge_primal_op.is_finite() and broken_hodge_dual_op.is_finite() );
            assert( broken_hodge_primal_op.getdimin()  == broken_mass_primal_op.getdimin() );
            assert( broken_hodge_primal_op.getdimout() == broken_mass_dual_op.getdimout()  );

            LOG << "sample vectors..." << nl;

            std::vector<FloatVector> trial_vectors;
            for( const auto& field : experiments_field[k] ) trial_vectors.push_back( Interpolation( M, n, k, r, field ) );

            for( int j = 0; j < number_of_random_trial_samples; j++ )
            {
                auto trial_vector = broken_mass_primal_op.createinputvector();
                trial_vector.random();
                trial_vector.normalize();
                assert( trial_vector.is_finite() );
                trial_vectors.push_back( trial_vector );
            }                
            


            /* We compare the following quantities      */
            /* 1. Mass squared u                        */
            /* 2. Mass squared of hodge_u               */
            /* 3. scalar integral of u vee u            */
            /* 4. volume integral of hodge_u wedge u    */
            /* 5. volume integral of u wedge hodge_u    */
            
            LOG << "testing..." << nl;
    
            errors_mass_dual [k][ l-l_min ][ r-r_min ] = -0.;
            errors_vee       [k][ l-l_min ][ r-r_min ] = -0.;
            errors_wedgeA    [k][ l-l_min ][ r-r_min ] = -0.;
            errors_wedgeB    [k][ l-l_min ][ r-r_min ] = -0.;

            sampled_mass_errors[k][l-l_min][r-r_min]   = -0.;

            
            for( int i = 0; i < trial_vectors.size(); i++ )
            {

                const auto&     field = trial_vectors[i];
                
                const auto dual_field = broken_hodge_primal_op * field;
                
                // mass( u, u )
                const Float mass_primal = field * ( broken_mass_primal_op * field );
                
                // mass( star u, star u )
                const Float mass_dual   = dual_field * ( broken_mass_dual_op * dual_field   );
                
                SparseMatrix broken_vee_with_primal_op   = FEECBrokenVeeMatrix(   M, M.getinnerdimension(),   k, r,   k, r, field );

                SparseMatrix broken_wedge_with_primal_op = FEECBrokenWedgeMatrix( M, M.getinnerdimension(), n-k, r,   k, r, field );

                SparseMatrix broken_wedge_with_dual_op   = FEECBrokenWedgeMatrix( M, M.getinnerdimension(),   k, r, n-k, r, dual_field );

                const auto vee_scalar = broken_vee_with_primal_op * field;

                const auto wedge_volume1 = broken_wedge_with_primal_op * dual_field;

                const auto wedge_volume2 = broken_wedge_with_dual_op   *      field;
                
                /* TODO(martin): be explicit about the sign switch ... no effect in 3D but an effect in other dimensions */

                // intscalar( u vee u )
                const Float vee_mass  = scalar_integrals * vee_scalar;
                
                // intvolume( u wedge star u )
                const Float wedge_mass1 = volume_integrals * wedge_volume1;
                
                // intvolume( star u wedge u )
                const Float wedge_mass2 = volume_integrals * wedge_volume2 * sign_power( k * (n-k) );
                
                LOGPRINTF("%.5e %.5e %.5e %.5e %.5e\n", mass_primal, mass_dual, vee_mass, wedge_mass1, wedge_mass2 );

                const auto error_mass   = absolute( mass_dual   - mass_primal );
                const auto error_vee    = absolute( vee_mass    - mass_primal );
                const auto error_wedgeA = absolute( wedge_mass1 - mass_primal );
                const auto error_wedgeB = absolute( wedge_mass2 - mass_primal );
                
                assert( field.is_finite()      );
                assert( dual_field.is_finite() );

                Assert( std::isfinite(error_mass) and std::isfinite(error_vee) and std::isfinite(error_wedgeA) and std::isfinite(error_wedgeB) );
                Assert( volume_integrals.getdimension() == wedge_volume1.getdimension(), volume_integrals.getdimension(), wedge_volume1.getdimension() );
                Assert( volume_integrals.getdimension() == wedge_volume2.getdimension(), volume_integrals.getdimension(), wedge_volume2.getdimension() );
                
                errors_mass_dual[k][l-l_min][r-r_min] = maximum( errors_mass_dual[k][l-l_min][r-r_min], error_mass   );
                errors_vee      [k][l-l_min][r-r_min] = maximum( errors_vee      [k][l-l_min][r-r_min], error_vee    );
                errors_wedgeA   [k][l-l_min][r-r_min] = maximum( errors_wedgeA   [k][l-l_min][r-r_min], error_wedgeA );
                errors_wedgeB   [k][l-l_min][r-r_min] = maximum( errors_wedgeB   [k][l-l_min][r-r_min], error_wedgeB );
                
                assert( std::isfinite(errors_mass_dual[k][l-l_min][r-r_min]) );
                assert( std::isfinite(errors_vee      [k][l-l_min][r-r_min]) );
                assert( std::isfinite(errors_wedgeA   [k][l-l_min][r-r_min]) );
                assert( std::isfinite(errors_wedgeB   [k][l-l_min][r-r_min]) );

                assert( field.is_finite() );
                assert( dual_field.is_finite() );

                LOG << "random test samples..." << nl;

                /* We verify the equality of:       */
                /* 1. mass(v,u)                     */
                /* 2. mass( star v, hodge_u )       */
                /* 3. mass( v, star hodge_u )       */
                /* 4. volint( hodge_u wedge v )     */
                    
                for( int j = 0; j < number_of_random_test_samples; j++ )
                {
                    auto test_field = broken_mass_primal_op.createinputvector();
                    test_field.random();
                    test_field.normalize();

                    // mass( v, u )
                    Float prod_test_mass_primal     = test_field * ( broken_mass_primal_op * field );
                    
                    // mass( hodge v, hodge_u )
                    Float prod_hodge_test_mass_dual = ( broken_hodge_primal_op * test_field ) * ( broken_mass_dual_op * dual_field );
                    
                    // mass( v, star hodge_u )
                    Float prod_test_mass_hodge_dual = ( broken_mass_primal_op * test_field ) * ( broken_hodge_dual_op * dual_field ) * sign_power( k * (n-k) );
                    
                    // intvolume( star u wedge v )
                    Float intvol_dual_wedge_test    = volume_integrals * ( broken_wedge_with_dual_op * test_field ) * sign_power( k * (n-k) );

                    // EXTRA: check self-inverse property
                    Float prod_trial_hodge_hodge_test = ( broken_mass_primal_op * field ) * ( ( broken_hodge_dual_op * broken_hodge_primal_op ) * test_field ) * sign_power( k * (n-k) );
                    
                    Float error_sampled = 0.;
                    error_sampled = maximum( error_sampled, absolute( prod_test_mass_primal - prod_hodge_test_mass_dual   ) );
                    error_sampled = maximum( error_sampled, absolute( prod_test_mass_primal - prod_test_mass_hodge_dual   ) );
                    error_sampled = maximum( error_sampled, absolute( prod_test_mass_primal - intvol_dual_wedge_test      ) );
                    error_sampled = maximum( error_sampled, absolute( prod_test_mass_primal - prod_trial_hodge_hodge_test ) );
                
                    LOGPRINTF("%d -- % .5e % .5e % .5e % .5e % .5e\n", j, prod_test_mass_primal, prod_hodge_test_mass_dual, prod_test_mass_hodge_dual, intvol_dual_wedge_test, prod_trial_hodge_hodge_test );
                    
                    assert( std::isfinite(prod_test_mass_primal) );
                    assert( std::isfinite(prod_hodge_test_mass_dual) );
                    assert( std::isfinite(prod_test_mass_hodge_dual) );
                    assert( std::isfinite(intvol_dual_wedge_test) );
                    assert( std::isfinite(prod_trial_hodge_hodge_test) );
                    assert( std::isfinite(error_sampled) );

                    sampled_mass_errors[k][l-l_min][r-r_min] = maximum( sampled_mass_errors[k][l-l_min][r-r_min], error_sampled );
                    
                    assert( std::isfinite(sampled_mass_errors[k][l-l_min][r-r_min]) );
                }

                /* Testing the alternative Hodge star matrix */ 

                auto alternative_dual_field = alternative_hodge_star * field;

                if(false)
                for( int j = 0; j < number_of_random_test_samples; j++ )
                {
                    auto test_dual_field = broken_mass_dual_op.createinputvector();
                    test_dual_field.random();
                    test_dual_field.normalize();

                    Float prod1 = test_dual_field * ( broken_mass_dual_op * alternative_dual_field );
                    Float prod2 = test_dual_field * ( broken_mass_dual_op *             dual_field );

                    LOGPRINTF("%d -- % .5e % .5e\n", j, prod1, prod2 );
                }
                    
            } // i
            
        } // k, r






        
        
        
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
    std::vector<ConvergenceTable> contables_wedgeA(n + 1);
    std::vector<ConvergenceTable> contables_wedgeB(n + 1);

    ConvergenceTable contable_sampled[n+1];
    
    // Names 

    for (int k = 0; k <= n; k++) {
        contables_mass[k].table_name   = "Rounding errors D2K_mass"   + std::to_string(k);
        contables_vee[k].table_name    = "Rounding errors D2K_vee"    + std::to_string(k);
        contables_wedgeA[k].table_name = "Rounding errors D2K_wedgeA" + std::to_string(k);
        contables_wedgeB[k].table_name = "Rounding errors D2K_wedgeB" + std::to_string(k);
    }

    contable_sampled[0].table_name = "Sampled errors scalar";
    contable_sampled[1].table_name = "Sampled errors vector";
    contable_sampled[2].table_name = "Sampled errors volume";

    // column names
    for (int k = 0; k <= n; k++) {
        for (int r = r_min; r <= r_max; r++) {
            contables_mass[k]   << ("R" + std::to_string(r));
            contables_vee[k]    << ("R" + std::to_string(r));
            contables_wedgeA[k] << ("R" + std::to_string(r));
            contables_wedgeB[k] << ("R" + std::to_string(r));
            contable_sampled[k] << ("R" + std::to_string(r));
        }
        contables_mass[k]   << nl;
        contables_vee[k]    << nl;
        contables_wedgeA[k] << nl;
        contables_wedgeB[k] << nl;
        contable_sampled[k] << nl; 
    }


    LOG << "Fill in data" << nl;
    
    for (int k = 0; k <= n; k++)
    for (int l = l_min; l <= l_max; l++) 
    {
        for (int r = r_min; r <= r_max; r++) {
            contables_mass[k]   << errors_mass_dual[k][l - l_min][r - r_min];
            contables_vee[k]    << errors_vee      [k][l - l_min][r - r_min];
            contables_wedgeA[k] << errors_wedgeA   [k][l - l_min][r - r_min];
            contables_wedgeB[k] << errors_wedgeB   [k][l - l_min][r - r_min];
            
            contable_sampled[k] << sampled_mass_errors[k][l-l_min][r-r_min];
        }
        contables_mass[k]   << nl;
        contables_vee[k]    << nl;
        contables_wedgeA[k] << nl;
        contables_wedgeB[k] << nl;
        
        contable_sampled[k] << nl; 
    }

    for (int k = 0; k <= n; k++) {
        contables_mass[k].lg();
        contables_vee[k].lg();
        contables_wedgeA[k].lg();
        contables_wedgeB[k].lg();
        
        contable_sampled[k].lg(); 

        LOG << "                   " << nl;
    }

    
    
    LOG << "Check that differences (experiments) are below: " << desired_closeness << nl;
    
    for( int l = l_min; l <= l_max; l++ ) 
    for( int r = r_min; r <= r_max; r++ ) 
    for( int k =     0; k <=     n; k++ ) 
    {
        Assert( errors_mass_dual[k][l-l_min][r-r_min] < desired_closeness, errors_mass_dual[k][l-l_min][r-r_min], desired_closeness, k, r, l );
        Assert( errors_vee      [k][l-l_min][r-r_min] < desired_closeness, errors_vee      [k][l-l_min][r-r_min], desired_closeness, k, r, l );
        Assert( errors_wedgeA   [k][l-l_min][r-r_min] < desired_closeness, errors_wedgeA   [k][l-l_min][r-r_min], desired_closeness, k, r, l );
        Assert( errors_wedgeB   [k][l-l_min][r-r_min] < desired_closeness, errors_wedgeB   [k][l-l_min][r-r_min], desired_closeness, k, r, l );
    }
    
    const Float threshold = 0.1;

    LOG << "Check that differences (sampled) are below: " << threshold << nl;
    
    for( int l = l_min; l <= l_max; l++ ) 
    for( int r = r_min; r <= r_max; r++ ) 
    for( int k =     0; k <=     n; k++ ) 
    {
        Assert( sampled_mass_errors[k][l-l_min][r-r_min] < threshold, sampled_mass_errors[k][l-l_min][r-r_min], threshold, k, r, l );
    }
    
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
