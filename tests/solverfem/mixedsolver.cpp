
#include <cmath>

#include <functional>

#include "../../base/include.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/math.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.inclwhitney.hpp"
#include "../../fem/utilities.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test: 2D Maxwell System" << nl;
    
    if(true){

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.getCoordinates().scale( 1.1 );
        
        M.check();
        
        M.automatic_dirichlet_flags();

        
        
        // const Float xfeq = 2.;
        // const Float yfeq = 2.;
        
        std::function<FloatVector(const FloatVector&)> experiment_sol = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 2 );
                return FloatVector({ 
                    bumpfunction(vec[0])*bumpfunction(vec[1]) //std::sin( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
                    ,
                    bumpfunction(vec[0])*bumpfunction(vec[1]) //std::sin( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
                });
            };
        

        std::function<FloatVector(const FloatVector&)> experiment_ndiv = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 2 );
                return FloatVector( { 
                    - bumpfunction_dev(vec[0]) * bumpfunction(vec[1]) - bumpfunction(vec[0]) * bumpfunction_dev(vec[1])
                });
            };
        

        std::function<FloatVector(const FloatVector&)> experiment_curl = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 2 );
                return FloatVector( { // - partial_y + partial_x
                    - bumpfunction(vec[0]) * bumpfunction_dev(vec[1]) + bumpfunction_dev(vec[0]) * bumpfunction(vec[1])
                });
            };
        

        std::function<FloatVector(const FloatVector&)> experiment_rhs = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 2 );
                return FloatVector({
                    -
                    bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1])
                    -
                    bumpfunction(vec[0])        * bumpfunction_devdev(vec[1])
                    ,
                    -
                    bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1])
                    -
                    bumpfunction(vec[0])        * bumpfunction_devdev(vec[1])
                    });
            };

            
            
            
        

        ConvergenceTable contable_sigma("Error: Sigma");
        ConvergenceTable contable_u    ("Error: u");
        ConvergenceTable contable_du   ("Error: du");
        ConvergenceTable contable_iter ("Iteration percentage");
        ConvergenceTable contable_time ("Time in milliseconds");
        ConvergenceTable contable_res  ("Residual");
        
        contable_sigma.print_rowwise_instead_of_columnwise = true;
        contable_u.print_rowwise_instead_of_columnwise     = true;
        contable_du.print_rowwise_instead_of_columnwise    = true;
        contable_iter.print_rowwise_instead_of_columnwise  = true;
        contable_time.print_rowwise_instead_of_columnwise  = true;
        contable_res.print_rowwise_instead_of_columnwise   = true;
        
        contable_sigma.display_convergence_rates = false;
        contable_u.display_convergence_rates     = false;
        contable_du.display_convergence_rates    = false;
        contable_iter.display_convergence_rates  = false;
        contable_time.display_convergence_rates  = false;
        contable_res.display_convergence_rates   = false;
        
        bool do_crmcsr       = true;
        bool do_crmcsr_dia   = false; // TODO(martin): fix bug  
        bool do_crmcsr_tb    = false; // TODO(martin): convergence deteriorates 
        bool do_crmcsr_ssor  = true;
        bool do_crmcpp       = true; 
        bool do_hersooblock  = true; 
        bool do_minresblock  = false; // comparable to Herzog-Soodhalter 
        bool do_systemherzog = true;
        bool do_sparseherzog = true;
        
        if( do_crmcsr )       { contable_sigma << "CRMcsr"; contable_u << "CRMcsr"; contable_du << "CRMcsr"; contable_iter << "CRMcsr"; contable_time << "CRMcsr"; contable_res << "CRMcsr"; } 
        if( do_crmcsr_dia )   { contable_sigma << "CRMcsrD"; contable_u << "CRMcsrD"; contable_du << "CRMcsrD"; contable_iter << "CRMcsrD"; contable_time << "CRMcsrD"; contable_res << "CRMcsrD"; } 
        if( do_crmcsr_tb )    { contable_sigma << "CRMcsrT"; contable_u << "CRMcsrT"; contable_du << "CRMcsrT"; contable_iter << "CRMcsrT"; contable_time << "CRMcsrT"; contable_res << "CRMcsrT"; } 
        if( do_crmcsr_ssor )  { contable_sigma << "CRMcsrS"; contable_u << "CRMcsrS"; contable_du << "CRMcsrS"; contable_iter << "CRMcsrS"; contable_time << "CRMcsrS"; contable_res << "CRMcsrS"; } 
        if( do_crmcpp )       { contable_sigma << "CRMcpp"; contable_u << "CRMcpp"; contable_du << "CRMcpp"; contable_iter << "CRMcpp"; contable_time << "CRMcpp"; contable_res << "CRMcpp"; } 
        if( do_hersooblock )  { contable_sigma << "Herzog"; contable_u << "Herzog"; contable_du << "Herzog"; contable_iter << "Herzog"; contable_time << "Herzog"; contable_res << "Herzog"; } 
        if( do_minresblock )  { contable_sigma << "Minres"; contable_u << "Minres"; contable_du << "Minres"; contable_iter << "Minres"; contable_time << "Minres"; contable_res << "Minres"; } 
        if( do_systemherzog ) { contable_sigma << "SysHerzog"; contable_u << "SysHerzog"; contable_du << "SysHerzog"; contable_iter << "SysHerzog"; contable_time << "SysHerzog"; contable_res << "SysHerzog"; } 
        if( do_sparseherzog ) { contable_sigma << "SpaHerzog"; contable_u << "SpaHerzog"; contable_du << "SpaHerzog"; contable_iter << "SpaHerzog"; contable_time << "SpaHerzog"; contable_res << "SpaHerzog"; } 
        
        { contable_sigma << nl; contable_u << nl; contable_du << nl; contable_iter << nl; contable_time << nl; contable_res << nl; } 
        

        const int min_l = 0; 
        
        const int max_l = 5;
        
        
        assert( 0 <= min_l and min_l <= max_l );
        
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ )
        {
            
            LOG << "Level: " << min_l << " <= " << l << " <= " << max_l << nl;
            LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            const int r = 1;

            {
                
                LOG << "... assemble partial matrices" << nl;
        
                SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r   );
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
                SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

                SparseMatrix vector_elevationmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r-1, 1);
                SparseMatrix vector_elevationmatrix_t = vector_elevationmatrix.getTranspose();

                SparseMatrix scalar_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r   );
                SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

                SparseMatrix vector_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r   );
                SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );
                SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();


                LOG << "... full matrices" << nl;
        
                auto mass = vector_incmatrix_t * vector_massmatrix * vector_incmatrix;

                auto mat_A  = scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix;
                mat_A.sortandcompressentries();
                
                LOG << vector_elevationmatrix.getdimin() << space << vector_elevationmatrix.getdimout() << nl;

                auto mat_Bt = scalar_incmatrix_t & scalar_diffmatrix_t & vector_elevationmatrix_t & vector_massmatrix & vector_incmatrix; // upper right
                mat_Bt.sortandcompressentries();
                
                auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower left
                mat_B.sortandcompressentries();
                
                auto mat_C  = vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix;
                mat_C.sortandcompressentries();
                
                auto A  = MatrixCSR( mat_A  );
                auto Bt = MatrixCSR( mat_Bt );
                auto B  = MatrixCSR( mat_B  );
                auto C  = MatrixCSR( mat_C  );
                
                auto PA = MatrixCSR( scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix )//IdentityOperator(A.getdimin())
                            + MatrixCSR( scalar_incmatrix_t & scalar_diffmatrix_t & vector_elevationmatrix_t & vector_massmatrix & vector_elevationmatrix & scalar_diffmatrix & scalar_incmatrix );
                auto PC = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )//IdentityOperator(C.getdimin())
                            + C;
                            
                auto SystemMatrix = C + B * inv(A,100*machine_epsilon) * Bt;
                
                {

                    LOG << "... interpolate explicit solution and rhs" << nl;
                    
                    FloatVector interpol_ndiv = Interpolation( M, M.getinnerdimension(), 0, r, experiment_ndiv  );
                    FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 1, r,   experiment_sol  );
                    FloatVector interpol_curl = Interpolation( M, M.getinnerdimension(), 2, r-1, experiment_curl );
                    FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 1, r,   experiment_rhs  );
                    FloatVector rhs = vector_incmatrix_t * ( vector_massmatrix * interpol_rhs );
                    
                    FloatVector sol( vector_incmatrix.getdimin(), 0. );

                    assert( do_crmcsr or do_crmcpp or do_hersooblock or do_minresblock or do_systemherzog or do_sparseherzog );
                    
                    for( int k = 0; k <= 8; k++ )
                    {

                        sol.zero();
                        
                        Float runtime = -1.;
                        int iteration_count = -1;

                        if( k == 0 and not do_crmcsr       ) continue;
                        if( k == 1 and not do_crmcsr_dia   ) continue;
                        if( k == 2 and not do_crmcsr_tb    ) continue;
                        if( k == 3 and not do_crmcsr_ssor  ) continue;
                        if( k == 4 and not do_crmcpp       ) continue;
                        if( k == 5 and not do_hersooblock  ) continue;
                        if( k == 6 and not do_minresblock  ) continue;
                        if( k == 7 and not do_systemherzog ) continue;
                        if( k == 8 and not do_sparseherzog ) continue;
                        
                        if( k == 0 and do_crmcsr )
                        {
                            sol.zero();
                            
                            FloatVector res = rhs;

                            timestamp start = timestampnow();

                            iteration_count = 
                            HodgeConjugateResidualSolverCSR(
                                B.getdimout(), 
                                A.getdimout(), 
                                sol.raw(), 
                                rhs.raw(), 
                                A.getA(),    A.getC(),    A.getV(), 
                                B.getA(),    B.getC(),    B.getV(), 
                                Bt.getA(),   Bt.getC(),   Bt.getV(), 
                                C.getA(),    C.getC(),    C.getV(),
                                res.raw(),
                                desired_precision,
                                100,
                                desired_precision,
                                -1
                            );
    
                            timestamp end = timestampnow();
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            runtime  = static_cast<Float>( end - start );
                        }


                        if( k == 1 and do_crmcsr_dia )
                        {
                            sol.zero();
                            
                            FloatVector res = rhs;

                            timestamp start = timestampnow();

                            iteration_count = 
                            HodgeConjugateResidualSolverCSR_diagonal(
                                B.getdimout(), 
                                A.getdimout(), 
                                sol.raw(), 
                                rhs.raw(), 
                                A.getA(),    A.getC(),    A.getV(), 
                                B.getA(),    B.getC(),    B.getV(), 
                                Bt.getA(),   Bt.getC(),   Bt.getV(), 
                                C.getA(),    C.getC(),    C.getV(),
                                res.raw(),
                                desired_precision,
                                100,
                                desired_precision,
                                -1
                            );
    
                            timestamp end = timestampnow();
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            runtime  = static_cast<Float>( end - start );
                        }


                        if( k == 2 and do_crmcsr_tb )
                        {
                            sol.zero();
                            
                            FloatVector res = rhs;

                            timestamp start = timestampnow();

                            iteration_count = 
                            HodgeConjugateResidualSolverCSR_textbook(
                                B.getdimout(), 
                                A.getdimout(), 
                                sol.raw(), 
                                rhs.raw(), 
                                A.getA(),    A.getC(),    A.getV(), 
                                B.getA(),    B.getC(),    B.getV(), 
                                Bt.getA(),   Bt.getC(),   Bt.getV(), 
                                C.getA(),    C.getC(),    C.getV(),
                                res.raw(),
                                desired_precision,
                                100,
                                desired_precision,
                                -1
                            );
    
                            timestamp end = timestampnow();
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            runtime  = static_cast<Float>( end - start );
                        }


                        if( k == 3 and do_crmcsr_ssor )
                        {
                            sol.zero();
                            
                            FloatVector res = rhs;

                            timestamp start = timestampnow();

                            iteration_count = 
                            HodgeConjugateResidualSolverCSR_SSOR(
                                B.getdimout(), 
                                A.getdimout(), 
                                sol.raw(), 
                                rhs.raw(), 
                                A.getA(),    A.getC(),    A.getV(), 
                                B.getA(),    B.getC(),    B.getV(), 
                                Bt.getA(),   Bt.getC(),   Bt.getV(), 
                                C.getA(),    C.getC(),    C.getV(),
                                res.raw(),
                                desired_precision,
                                100,
                                desired_precision,
                                -1
                            );
    
                            timestamp end = timestampnow();
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            runtime  = static_cast<Float>( end - start );
                        }


                        if( k == 4 and do_crmcpp )
                        {
                            sol.zero();
                            
                            ConjugateResidualMethod solver( SystemMatrix );
                            solver.print_modulo        = 100;
                            solver.max_iteration_count = 1 * sol.getdimension();

                            timestamp start = timestampnow();
                            solver.solve( sol, rhs );
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            runtime  = static_cast<Float>( end - start );

                            iteration_count = solver.recent_iteration_count;
                        }


                        if( k == 5 and do_hersooblock )
                        {
                            sol.zero();
                            
                            auto X = Block2x2Operator( A.getdimout() + B.getdimout(), A.getdimin() + Bt.getdimin(), -A, Bt, B, C );

                            FloatVector sol_whole( A.getdimin()  + Bt.getdimin(),  0. );
                            FloatVector rhs_whole( A.getdimout() +  B.getdimout(), 0. );
                            
                            sol_whole.zero();
                            rhs_whole.setslice( 0, A.getdimout(), 0. );
                            rhs_whole.setslice( A.getdimout(), rhs );
                            
                            HerzogSoodhalterMethod solver( X );
                            solver.print_modulo        = 100;
                            solver.max_iteration_count = 1 * sol_whole.getdimension();

                            timestamp start = timestampnow();
                            solver.solve( sol_whole, rhs_whole );
                            timestamp end = timestampnow();

                            sol = sol_whole.getslice( A.getdimout(), B.getdimout() );
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            runtime  = static_cast<Float>( end - start );

                            iteration_count = solver.recent_iteration_count;
                        }


                        if( k == 6 and do_minresblock )
                        {   
                            sol.zero();
                            
                            auto X = Block2x2Operator( A.getdimout() + B.getdimout(), A.getdimin() + Bt.getdimin(), -A, Bt, B, C );

                            FloatVector sol_whole( A.getdimin()  + Bt.getdimin(),  0. );
                            FloatVector rhs_whole( A.getdimout() +  B.getdimout(), 0. );
                            
                            sol_whole.zero();
                            rhs_whole.setslice( 0, A.getdimout(), 0. );
                            rhs_whole.setslice( A.getdimout(), rhs );
                            
                            MinimumResidualMethod solver( X );
                            solver.print_modulo        = 100;
                            solver.max_iteration_count = 1 * sol_whole.getdimension();

                            timestamp start = timestampnow();
                            solver.solve( sol_whole, rhs_whole );
                            timestamp end = timestampnow();

                            sol = sol_whole.getslice( A.getdimout(), B.getdimout() );
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            runtime  = static_cast<Float>( end - start );

                            iteration_count = solver.recent_iteration_count;
                        }


                        if( k == 7 and do_systemherzog )
                        {   
                            sol.zero();
                            
                            // const auto PAinv = IdentityOperator(A.getdimin());
                            // const auto PCinv = IdentityOperator(C.getdimin());
                            const auto PAinv = inv(PA,desired_precision,0);
                            const auto PCinv = inv(PC,desired_precision,0);

                            FloatVector  x_A( A.getdimin(),  0. ); 
                            FloatVector& x_C = sol;
                            
                            const FloatVector  b_A( A.getdimin(),  0. ); 
                            const FloatVector& b_C = rhs; 
                            
                            timestamp start = timestampnow();
                            iteration_count = 
                            BlockHerzogSoodhalterMethod( 
                                x_A, 
                                x_C, 
                                b_A, 
                                b_C, 
                                -A, Bt, B, C, 
                                desired_precision,
                                1,
                                PAinv, PCinv
                            );
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            runtime  = static_cast<Float>( end - start );
                        }


                        if( k == 8 and do_sparseherzog )
                        {   
                            sol.zero();
                            
                            FloatVector  x_A( A.getdimin(),  0. ); 
                            FloatVector& x_C = sol;
                            
                            const FloatVector  b_A( A.getdimin(),  0. ); 
                            const FloatVector& b_C = rhs; 
                            
                            timestamp start = timestampnow();
                            iteration_count = 
                            HodgeHerzogSoodhalterMethod( 
                                x_A.getdimension(), 
                                x_C.getdimension(), 
                                x_A.raw(), 
                                x_C.raw(), 
                                b_A.raw(),  
                                b_C.raw(), 
                                A.getA(),    A.getC(),    A.getV(), 
                                B.getA(),    B.getC(),    B.getV(), 
                                Bt.getA(),   Bt.getC(),   Bt.getV(), 
                                C.getA(),    C.getC(),    C.getV(),
                                desired_precision,
                                1,
                                // nullptr, nullptr, nullptr,
                                // nullptr, nullptr, nullptr,
                                PA.getA(), PA.getC(), PA.getV(),
                                PC.getA(), PC.getC(), PC.getV(),
                                desired_precision, 
                                0
                            );
                            
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            runtime  = static_cast<Float>( end - start );
                        }

                        assert( runtime >= 0. and iteration_count >= 0 );

                        assert( sol.is_finite() );

                        auto ndiv = inv(A,desired_precision) * Bt * sol;
                        
                        auto curl = vector_diffmatrix * vector_incmatrix * sol;
                        
                        LOG << "... compute error and residual:" << k << nl;

                        auto errornorm_aux_ndiv = interpol_ndiv - scalar_incmatrix * ndiv;
                        auto errornorm_aux_sol  = interpol_sol  - vector_incmatrix *  sol;
                        auto errornorm_aux_curl = interpol_curl -                    curl;

                        Float errornorm_ndiv = std::sqrt( errornorm_aux_ndiv * ( scalar_massmatrix * errornorm_aux_ndiv ) );
                        Float errornorm_sol  = std::sqrt( errornorm_aux_sol  * ( vector_massmatrix * errornorm_aux_sol  ) );
                        Float errornorm_curl = std::sqrt( errornorm_aux_curl * ( volume_massmatrix * errornorm_aux_curl ) );
                        Float residualnorm   = ( rhs - SystemMatrix * sol ).norm();
                        
                        contable_sigma << errornorm_ndiv;
                        contable_u     << errornorm_sol;
                        contable_du    << errornorm_curl;
                        contable_res   << residualnorm;
                        contable_iter  << iteration_count / static_cast<Float>( A.getdimin() + C.getdimin() );
                        contable_time  << runtime;
                        
                    }

                    contable_sigma << nl;
                    contable_u     << nl;
                    contable_du    << nl;
                    contable_res   << nl;
                    contable_iter  << nl;
                    contable_time  << nl;
                
                        
                    contable_sigma.lg();
                    contable_u    .lg();
                    contable_du   .lg();
                    contable_res  .lg();
                    contable_iter .lg();
                    contable_time .lg();
                    


                }
                
            }

            if( l != max_l ) { 
                LOG << "Refinement..." << nl; 
                M.uniformrefinement();
                M.shake_interior_vertices();
            }
    
        } 
        
        contable_sigma.lg();
        contable_u    .lg();
        contable_du   .lg();
        contable_res  .lg();
        contable_iter .lg();
        contable_time .lg();
    
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}




