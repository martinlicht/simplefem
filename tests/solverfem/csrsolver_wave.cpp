

/**/

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
    LOG << "Unit Test: 2D Maxwell System" << endl;
    
    // LOG << std::setprecision(10);

    LOG << "Initial mesh..." << endl;
    
    MeshSimplicial2D M = StandardSquare2D();
    
    M.getcoordinates().scale( 1.1 );
    
    M.check();
    
    M.automatic_dirichlet_flags();

    
    
    const Float xfeq = 2.;
    const Float yfeq = 2.;
    
    std::function<FloatVector(const FloatVector&)> experiment_sol = 
        [=](const FloatVector& vec) -> FloatVector{
            Assert( vec.getdimension() == 2 );
            return FloatVector({ 
                bumpfunction(vec[0])*bumpfunction(vec[1]) //std::sin( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
                ,
                bumpfunction(vec[0])*bumpfunction(vec[1]) //std::sin( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
            });
        };
    
    std::function<FloatVector(const FloatVector&)> experiment_curl = 
        [=](const FloatVector& vec) -> FloatVector{
            Assert( vec.getdimension() == 2 );
            return FloatVector( { // - partial_y + partial_x
                - bumpfunction(vec[0]) * bumpfunction_dev(vec[1]) + bumpfunction_dev(vec[0]) * bumpfunction(vec[1])
            });
        };

    std::function<FloatVector(const FloatVector&)> experiment_rhs = 
        [=](const FloatVector& vec) -> FloatVector{
            Assert( vec.getdimension() == 2 );
            return FloatVector({
                bumpfunction(vec[0])*bumpfunction(vec[1])
                -
                bumpfunction(vec[0]) * bumpfunction_devdev(vec[1])
                +
                bumpfunction_dev(vec[0]) * bumpfunction_dev(vec[1])
                ,
                bumpfunction(vec[0])*bumpfunction(vec[1])
                -
                bumpfunction_devdev(vec[0]) * bumpfunction(vec[1])
                +
                bumpfunction_dev(vec[0]) * bumpfunction_dev(vec[1])
                });
        };

        
        
        
    

    ConvergenceTable contable_u    ("Error: u");
    ConvergenceTable contable_du   ("Error: du");
    ConvergenceTable contable_iter ("Iteration percentage");
    ConvergenceTable contable_time ("Time in milliseconds");
    ConvergenceTable contable_res  ("Residual");
    
    contable_u.print_rowwise_instead_of_columnwise     = true;
    contable_du.print_rowwise_instead_of_columnwise    = true;
    contable_iter.print_rowwise_instead_of_columnwise  = true;
    contable_time.print_rowwise_instead_of_columnwise  = true;
    contable_res.print_rowwise_instead_of_columnwise   = true;
    
    contable_u.display_convergence_rates     = false;
    contable_du.display_convergence_rates    = false;
    contable_iter.display_convergence_rates  = false;
    contable_time.display_convergence_rates  = false;
    contable_res.display_convergence_rates   = false;
    

    bool do_cgmpp      = true;
    bool do_crmpp_expl = true;
    bool do_crmpp_robt = true;
    bool do_crmpp_fast = true;
    bool do_minres     = true;
    bool do_herzog     = true;
    //
    bool do_cgm_csr                = true;
    bool do_crm_csr                = true;
    bool do_crm_csrtextbook        = true;
    bool do_minres_csr             = true;
    bool do_whatever_csr           = false;
    bool do_cgm_diagonal_csr       = true;
    bool do_cgm_ssor_csr           = true;
    bool do_chebyshev_diagonal_csr = false;

    if( do_cgmpp      ) contable_u << "CGM++"      ;
    if( do_crmpp_expl ) contable_u << "CRM++(expl)";
    if( do_crmpp_robt ) contable_u << "CRM++(robt)";
    if( do_crmpp_fast ) contable_u << "CRM++(fast)";
    if( do_minres     ) contable_u << "MINRES"     ;
    if( do_herzog     ) contable_u << "HERZOG"     ;
    //
    if( do_cgm_csr )                contable_u << "CGMcsr"       ;
    if( do_crm_csr )                contable_u << "CRMcsr"       ;
    if( do_crm_csrtextbook )        contable_u << "CRMcsr_tb"    ;
    if( do_minres_csr )             contable_u << "MINREScsr"    ;
    if( do_whatever_csr )           contable_u << "WHATEVER"     ;
    if( do_cgm_diagonal_csr )       contable_u << "CGMcsr_diag"  ;
    if( do_cgm_ssor_csr )           contable_u << "CGMcsr_ssor"  ;
    if( do_chebyshev_diagonal_csr ) contable_u << "Chebyshev_csr";
    
    if( do_cgmpp      ) contable_du << "CGM++"      ;
    if( do_crmpp_expl ) contable_du << "CRM++(expl)";
    if( do_crmpp_robt ) contable_du << "CRM++(robt)";
    if( do_crmpp_fast ) contable_du << "CRM++(fast)";
    if( do_minres     ) contable_du << "MINRES"     ;
    if( do_herzog     ) contable_du << "HERZOG"     ;
    //
    if( do_cgm_csr )                contable_du << "CGMcsr"       ;
    if( do_crm_csr )                contable_du << "CRMcsr"       ;
    if( do_crm_csrtextbook )        contable_du << "CRMcsr_tb"    ;
    if( do_minres_csr )             contable_du << "MINREScsr"    ;
    if( do_whatever_csr )           contable_du << "WHATEVER"     ;
    if( do_cgm_diagonal_csr )       contable_du << "CGMcsr_diag"  ;
    if( do_cgm_ssor_csr )           contable_du << "CGMcsr_ssor"  ;
    if( do_chebyshev_diagonal_csr ) contable_du << "Chebyshev_csr";
    
    if( do_cgmpp      ) contable_iter << "CGM++"      ;
    if( do_crmpp_expl ) contable_iter << "CRM++(expl)";
    if( do_crmpp_robt ) contable_iter << "CRM++(robt)";
    if( do_crmpp_fast ) contable_iter << "CRM++(fast)";
    if( do_minres     ) contable_iter << "MINRES"     ;
    if( do_herzog     ) contable_iter << "HERZOG"     ;
    //
    if( do_cgm_csr )                contable_iter << "CGMcsr"       ;
    if( do_crm_csr )                contable_iter << "CRMcsr"       ;
    if( do_crm_csrtextbook )        contable_iter << "CRMcsr_tb"    ;
    if( do_minres_csr )             contable_iter << "MINREScsr"    ;
    if( do_whatever_csr )           contable_iter << "WHATEVER"     ;
    if( do_cgm_diagonal_csr )       contable_iter << "CGMcsr_diag"  ;
    if( do_cgm_ssor_csr )           contable_iter << "CGMcsr_ssor"  ;
    if( do_chebyshev_diagonal_csr ) contable_iter << "Chebyshev_csr";

    if( do_cgmpp      ) contable_time << "CGM++"      ;
    if( do_crmpp_expl ) contable_time << "CRM++(expl)";
    if( do_crmpp_robt ) contable_time << "CRM++(robt)";
    if( do_crmpp_fast ) contable_time << "CRM++(fast)";
    if( do_minres     ) contable_time << "MINRES"     ;
    if( do_herzog     ) contable_time << "HERZOG"     ;
    //
    if( do_cgm_csr )                contable_time << "CGMcsr"       ;
    if( do_crm_csr )                contable_time << "CRMcsr"       ;
    if( do_crm_csrtextbook )        contable_time << "CRMcsr_tb"    ;
    if( do_minres_csr )             contable_time << "MINREScsr"    ;
    if( do_whatever_csr )           contable_time << "WHATEVER"     ;
    if( do_cgm_diagonal_csr )       contable_time << "CGMcsr_diag"  ;
    if( do_cgm_ssor_csr )           contable_time << "CGMcsr_ssor"  ;
    if( do_chebyshev_diagonal_csr ) contable_time << "Chebyshev_csr";

    if( do_cgmpp      ) contable_res << "CGM++"      ;
    if( do_crmpp_expl ) contable_res << "CRM++(expl)";
    if( do_crmpp_robt ) contable_res << "CRM++(robt)";
    if( do_crmpp_fast ) contable_res << "CRM++(fast)";
    if( do_minres     ) contable_res << "MINRES"     ;
    if( do_herzog     ) contable_res << "HERZOG"     ;
    //
    if( do_cgm_csr )                contable_res << "CGMcsr"       ;
    if( do_crm_csr )                contable_res << "CRMcsr"       ;
    if( do_crm_csrtextbook )        contable_res << "CRMcsr_tb"    ;
    if( do_minres_csr )             contable_res << "MINREScsr"    ;
    if( do_whatever_csr )           contable_res << "WHATEVER"     ;
    if( do_cgm_diagonal_csr )       contable_res << "CGMcsr_diag"  ;
    if( do_cgm_ssor_csr )           contable_res << "CGMcsr_ssor"  ;
    if( do_chebyshev_diagonal_csr ) contable_res << "Chebyshev_csr";

    const int min_l = 0; 
    
    const int max_l = 6;
    
    const int min_r = 1; 
    
    const int max_r = 1;

    assert( 0 <= min_l and min_l <= max_l );
    assert( 0 <= min_r and min_r <= max_r );
    
    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ )
    {
        
        LOG << "Level: " << l << "/" << max_l << std::endl;
        LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        for( int r = min_r; r <= max_r; r++ )
        {
            
            LOG << "Polynomial degree: " << r << "/" << max_r << std::endl;
            
            LOG << "...assemble partial matrices" << endl;
    
            SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

            SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

            SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
            SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

            LOG << "...convert to CSR" << endl;
    
            auto MassMatrix = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix );
            auto DiffMatrix = MatrixCSR( vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix );
            auto SystemMatrix = MassMatrix + DiffMatrix; 
            
            {

                LOG << "...interpolate explicit solution and rhs" << endl;
                
                FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 1, r,   experiment_sol  );
                FloatVector interpol_curl = Interpolation( M, M.getinnerdimension(), 2, r-1, experiment_curl );
                FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 1, r,   experiment_rhs  );
                
                const FloatVector rhs = vector_incmatrix_t * ( vector_massmatrix * interpol_rhs );
                
                FloatVector sol_orignal( SystemMatrix.getdimin(), 0. );

                for( int k = 0; k <= 12; k++ )
                {

                    Float runtime;
                    int iteration_count;

                    if( k== 0 and not do_cgmpp      ) continue;
                    if( k== 1 and not do_crmpp_expl ) continue;
                    if( k== 2 and not do_crmpp_robt ) continue;
                    if( k== 3 and not do_crmpp_fast ) continue;
                    if( k== 4 and not do_minres     ) continue;
                    if( k== 5 and not do_herzog     ) continue;
                    //
                    if( k== 6 and not do_cgm_csr )                continue;
                    if( k== 7 and not do_crm_csr )                continue;
                    if( k== 8 and not do_crm_csrtextbook )        continue;
                    if( k== 9 and not do_minres_csr )             continue;
                    if( k==10 and not do_whatever_csr )           continue;
                    if( k==11 and not do_cgm_diagonal_csr )       continue;
                    if( k==12 and not do_cgm_ssor_csr )           continue;
                    if( k==13 and not do_chebyshev_diagonal_csr ) continue;

                    auto sol = sol_orignal;
                        

                    if( k== 0 and do_cgmpp )
                    {
                        LOG << "CGM C++" << endl;
                    
                        ConjugateGradientMethod Solver( SystemMatrix );
                        Solver.print_modulo        = 0;
                        Solver.threshold           = desired_precision;
                        Solver.max_iteration_count = 1 * sol.getdimension();
                        timestamp start = gettimestamp();
                        Solver.solve( sol, rhs );
                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        runtime  = static_cast<Float>( end - start );
                        iteration_count = Solver.recent_iteration_count;
                    }

                    if( k== 1 and do_crmpp_expl )
                    {
                        LOG << "CRM C++" << endl;
                    
                        ConjugateResidualMethod Solver( SystemMatrix );
                        Solver.print_modulo        = 0;
                        Solver.threshold           = desired_precision;
                        Solver.max_iteration_count = 1 * sol.getdimension();
                        timestamp start = gettimestamp();
                        Solver.solve_explicit( sol, rhs );
                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        runtime  = static_cast<Float>( end - start );
                        iteration_count = Solver.recent_iteration_count;
                    }

                    if( k== 2 and do_crmpp_robt )
                    {
                        LOG << "CRM C++" << endl;
                    
                        ConjugateResidualMethod Solver( SystemMatrix );
                        Solver.print_modulo        = 0;
                        Solver.threshold           = desired_precision;
                        Solver.max_iteration_count = 1 * sol.getdimension();
                        timestamp start = gettimestamp();
                        Solver.solve_robust( sol, rhs );
                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        runtime  = static_cast<Float>( end - start );
                        iteration_count = Solver.recent_iteration_count;
                    }

                    if( k== 3 and do_crmpp_fast )
                    {
                        LOG << "CRM C++" << endl;
                    
                        ConjugateResidualMethod Solver( SystemMatrix );
                        Solver.print_modulo        = 0;
                        Solver.threshold           = desired_precision;
                        Solver.max_iteration_count = 1 * sol.getdimension();
                        timestamp start = gettimestamp();
                        Solver.solve_fast( sol, rhs );
                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        runtime  = static_cast<Float>( end - start );
                        iteration_count = Solver.recent_iteration_count;
                    }

                    if( k== 4 and do_minres )
                    {
                        LOG << "MINRES C++" << endl;
                    
                        MinimumResidualMethod Solver( SystemMatrix );
                        Solver.print_modulo        = 0;
                        Solver.threshold           = desired_precision;
                        Solver.max_iteration_count = 1 * sol.getdimension();
                        timestamp start = gettimestamp();
                        Solver.solve( sol, rhs );
                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                        runtime  = static_cast<Float>( end - start );
                        iteration_count = Solver.recent_iteration_count;
                    }

                    if( k== 5 and do_herzog )
                    {
                        LOG << "HERZOG SOODHALTER C++" << endl;
                    
                        HerzogSoodhalterMethod Solver( SystemMatrix );
                        Solver.print_modulo        = 0;
                        Solver.threshold           = desired_precision;
                        Solver.max_iteration_count = 1 * sol.getdimension();
                        timestamp start = gettimestamp();
                        Solver.solve( sol, rhs );
                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                        runtime  = static_cast<Float>( end - start );
                        iteration_count = Solver.recent_iteration_count;
                    }
                    
                    if( k== 6 and do_cgm_csr )
                    {
                        LOG << "CGM - CSR Classic" << endl;
                    
                        FloatVector residual( rhs );
                        auto max_iteration_count = sol.getdimension();
                        timestamp start = gettimestamp();
                        auto recent_iteration_count = 
                        ConjugateGradientSolverCSR( 
                            sol.getdimension(), 
                            sol.raw(), 
                            rhs.raw(), 
                            SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
                            residual.raw(),
                            desired_precision,
                            0
                        );
                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        runtime  = static_cast<Float>( end - start );
                        iteration_count = recent_iteration_count;
                    }

                    if( k== 7 and do_crm_csr )
                    {
                        LOG << "CRM - CSR Classic" << endl;
                    
                        FloatVector residual( rhs );
                        auto max_iteration_count = sol.getdimension();
                        timestamp start = gettimestamp();
                        auto recent_iteration_count = 
                        ConjugateResidualSolverCSR( 
                            sol.getdimension(), 
                            sol.raw(), 
                            rhs.raw(), 
                            SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
                            residual.raw(),
                            desired_precision,
                            0
                        );

                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        runtime  = static_cast<Float>( end - start );
                        iteration_count = recent_iteration_count;
                    }

                    if( k== 8 and do_crm_csrtextbook )
                    {
                        LOG << "CRM - CSR Textbook" << endl;
                    
                        FloatVector residual( rhs );
                        auto max_iteration_count = sol.getdimension();
                        timestamp start = gettimestamp();
                        auto recent_iteration_count = 
                        ConjugateResidualSolverCSR_textbook( 
                            sol.getdimension(), 
                            sol.raw(), 
                            rhs.raw(), 
                            SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
                            residual.raw(),
                            desired_precision,
                            0
                        );

                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        runtime  = static_cast<Float>( end - start );
                        iteration_count = recent_iteration_count;
                    }

                    if( k== 9 and do_minres_csr )
                    {
                        LOG << "MINRES CSR" << endl;
                    
                        FloatVector residual( rhs );
                        auto max_iteration_count = sol.getdimension();
                        timestamp start = gettimestamp();
                        auto recent_iteration_count = 
                        MINRESCSR( 
                            sol.getdimension(), 
                            sol.raw(), 
                            rhs.raw(), 
                            SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
                            residual.raw(),
                            desired_precision,
                            0
                        );

                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        runtime  = static_cast<Float>( end - start );
                        iteration_count = recent_iteration_count;
                    }


                    if( k==10 and do_whatever_csr )
                    {
                        LOG << "WHATEVER CSR" << endl;
                    
                        FloatVector residual( rhs );
                        auto max_iteration_count = sol.getdimension();
                        timestamp start = gettimestamp();
                        auto recent_iteration_count = 
                        WHATEVER( 
                            sol.getdimension(), 
                            sol.raw(), 
                            rhs.raw(), 
                            SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
                            residual.raw(),
                            desired_precision,
                            0
                        );

                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        runtime  = static_cast<Float>( end - start );
                        iteration_count = recent_iteration_count;
                    }

                    if( k==11 and do_cgm_diagonal_csr )
                    {
                        LOG << "CGM - CSR Classic with diagonal preconditioning" << endl;
                        
                        auto precon = InverseDiagonalPreconditioner( SystemMatrix );

                        FloatVector residual( rhs );
                        auto max_iteration_count = sol.getdimension();
                        timestamp start = gettimestamp();
                        auto recent_iteration_count = 
                        ConjugateGradientSolverCSR_DiagonalPreconditioner( 
                            sol.getdimension(), 
                            sol.raw(), 
                            rhs.raw(), 
                            SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
                            residual.raw(),
                            desired_precision,
                            0,
                            precon.getdiagonal().raw()
                        );

                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        runtime  = static_cast<Float>( end - start );
                        iteration_count = recent_iteration_count;
                    }

                    if( k==12 and do_cgm_ssor_csr )
                    {
                        LOG << "CGM - CSR Classic with SSOR" << endl;
                        
                        auto diagonal = SystemMatrix.diagonal();

                        FloatVector residual( rhs );
                        auto max_iteration_count = sol.getdimension();
                        timestamp start = gettimestamp();
                        auto recent_iteration_count = 
                        ConjugateGradientSolverCSR_SSOR( 
                            sol.getdimension(), 
                            sol.raw(), 
                            rhs.raw(), 
                            SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
                            residual.raw(),
                            desired_precision,
                            0,
                            diagonal.raw(),
                            1.0
                        );

                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        runtime  = static_cast<Float>( end - start );
                        iteration_count = recent_iteration_count;
                    }

                    assert( sol.isfinite() );

                    LOG << "...compute error and residual:" << k << endl;

                    auto errornorm_aux_sol  = interpol_sol  - vector_incmatrix *  sol;
                    auto errornorm_aux_curl = interpol_curl - vector_diffmatrix * vector_incmatrix * sol;

                    Float errornorm_sol  = sqrt( errornorm_aux_sol  * ( vector_massmatrix * errornorm_aux_sol  ) );
                    Float errornorm_curl = sqrt( errornorm_aux_curl * ( volume_massmatrix * errornorm_aux_curl ) );
                    Float residualnorm   = ( rhs - SystemMatrix * sol ).norm();
                    
                    contable_u     << errornorm_sol;
                    contable_du    << errornorm_curl;
                    contable_res   << residualnorm;
                    contable_iter  << iteration_count / static_cast<Float>( SystemMatrix.getdimin() );
                    contable_time  << runtime;
                        
                }

                contable_u     << nl;
                contable_du    << nl;
                contable_res   << nl;
                contable_iter  << nl;
                contable_time  << nl;

                contable_u    .lg();
                contable_du   .lg();
                contable_res  .lg();
                contable_iter .lg();
                contable_time .lg();
                
            }
            
        }

        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }

    } 
    
    contable_u    .lg();
    contable_du   .lg();
    contable_res  .lg();
    contable_iter .lg();
    contable_time .lg();
    
    LOG << "Finished Unit Test" << endl;
    
    return 0;
}




