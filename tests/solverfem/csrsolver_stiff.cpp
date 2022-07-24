

/**/

#include <ostream>
// #include <fstream>
// #include <iomanip>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../operators/composedoperators.hpp"
// #include "../../operators/composed.hpp"
#include "../../dense/densematrix.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/sparsesolver.hpp"
// #include "../../solver/chebyshev.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: Compare numerical solvers CRM vs MINRES\n           for Solution of Dirichlet Problem" << endl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            M.automatic_dirichlet_flags();
            
            M.check_dirichlet_flags();
            
            LOG << "Prepare scalar fields for testing..." << endl;
            

            // std::function<FloatVector(const FloatVector&)> constant_one
            //     = [](const FloatVector& vec) -> FloatVector{
            //             assert( vec.getdimension() == 2 );
            //             return FloatVector({ 1. });
            //         };
            
            const Float xfeq = 1.;
            const Float yfeq = 1.;
            
            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        
                        xfeq*xfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] )
                        +
                        yfeq*yfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] )
                     });
                };
            

            
            
            

            LOG << "Solving Poisson Problem with Dirichlet boundary conditions" << endl;

            // ConvergenceTable contable_sol("L2 Error");
            ConvergenceTable contable_res("L2 Residual");
            ConvergenceTable contable_num("Iteration percentage");
            ConvergenceTable contable_sec("Time");

            // contable_sol.print_rowwise_instead_of_columnwise = true;
            contable_res.print_rowwise_instead_of_columnwise = true;
            contable_num.print_rowwise_instead_of_columnwise = true;
            contable_sec.print_rowwise_instead_of_columnwise = true;
            
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

            // if( do_cgmpp      ) contable_sol << "CGM++"      ;
            // if( do_crmpp_expl ) contable_sol << "CRM++(expl)";
            // if( do_crmpp_robt ) contable_sol << "CRM++(robt)";
            // if( do_crmpp_fast ) contable_sol << "CRM++(fast)";
            // if( do_minres     ) contable_sol << "MINRES"     ;
            // if( do_herzog     ) contable_sol << "HERZOG"     ;
            // //
            // if( do_cgm_csr )                contable_sol << "CGMcsr"       ;
            // if( do_crm_csr )                contable_sol << "CRMcsr"       ;
            // if( do_crm_csrtextbook )        contable_sol << "CRMcsr_tb"    ;
            // if( do_minres_csr )             contable_sol << "MINREScsr"    ;
            // if( do_whatever_csr )           contable_sol << "WHATEVER"     ;
            // if( do_cgm_diagonal_csr )       contable_sol << "CGMcsr_diag"  ;
            // if( do_cgm_ssor_csr )           contable_sol << "CGMcsr_ssor"  ;
            // if( do_chebyshev_diagonal_csr ) contable_sol << "Chebyshev_csr";
            
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

            if( do_cgmpp      ) contable_num << "CGM++"      ;
            if( do_crmpp_expl ) contable_num << "CRM++(expl)";
            if( do_crmpp_robt ) contable_num << "CRM++(robt)";
            if( do_crmpp_fast ) contable_num << "CRM++(fast)";
            if( do_minres     ) contable_num << "MINRES"     ;
            if( do_herzog     ) contable_num << "HERZOG"     ;
            //
            if( do_cgm_csr )                contable_num << "CGMcsr"       ;
            if( do_crm_csr )                contable_num << "CRMcsr"       ;
            if( do_crm_csrtextbook )        contable_num << "CRMcsr_tb"    ;
            if( do_minres_csr )             contable_num << "MINREScsr"    ;
            if( do_whatever_csr )           contable_num << "WHATEVER"     ;
            if( do_cgm_diagonal_csr )       contable_num << "CGMcsr_diag"  ;
            if( do_cgm_ssor_csr )           contable_num << "CGMcsr_ssor"  ;
            if( do_chebyshev_diagonal_csr ) contable_num << "Chebyshev_csr";

            if( do_cgmpp      ) contable_sec << "CGM++"      ;
            if( do_crmpp_expl ) contable_sec << "CRM++(expl)";
            if( do_crmpp_robt ) contable_sec << "CRM++(robt)";
            if( do_crmpp_fast ) contable_sec << "CRM++(fast)";
            if( do_minres     ) contable_sec << "MINRES"     ;
            if( do_herzog     ) contable_sec << "HERZOG"     ;
            //
            if( do_cgm_csr )                contable_sec << "CGMcsr"       ;
            if( do_crm_csr )                contable_sec << "CRMcsr"       ;
            if( do_crm_csrtextbook )        contable_sec << "CRMcsr_tb"    ;
            if( do_minres_csr )             contable_sec << "MINREScsr"    ;
            if( do_whatever_csr )           contable_sec << "WHATEVER"     ;
            if( do_cgm_diagonal_csr )       contable_sec << "CGMcsr_diag"  ;
            if( do_cgm_ssor_csr )           contable_sec << "CGMcsr_ssor"  ;
            if( do_chebyshev_diagonal_csr ) contable_sec << "Chebyshev_csr";
            

            const int min_l = 3;
            
            const int max_l = 7;

            assert( 0 <= min_l and min_l <= max_l );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << "/" << max_l << std::endl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                const int r = 1;
                
                {
                    
                    LOG << "...assemble scalar mass matrix" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    LOG << "...assemble vector mass matrix" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );

                    LOG << "...sort and compress mass matrices" << endl;
            
                    scalar_massmatrix.sortandcompressentries();
                    vector_massmatrix.sortandcompressentries();

                    LOG << "...assemble differential matrix and transpose" << endl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    LOG << "...assemble inclusion matrix and transpose" << endl;
            
                    SparseMatrix incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    LOG << "...compose stiffness and mass matrices" << endl;
            
                    const auto composed_stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
                    const auto composed_mass      = incmatrix_t * scalar_massmatrix * incmatrix;

                    auto opr  = diffmatrix & incmatrix;
                    auto opl  = opr.getTranspose(); 
                    auto stiffness_csr_prelim = opl & ( vector_massmatrix & opr );
                    
                    LOG << "...convert to CSR" << endl;
            
                    stiffness_csr_prelim.sortentries();
                    auto stiffness_csr = MatrixCSR( stiffness_csr_prelim );

                    const auto& stiffness = stiffness_csr;
                    const auto& mass      = composed_mass;
                    
                    LOG << "...matrices done" << endl;

                    {

                        FloatVector sol_original( M.count_simplices(0), 0. );
                        sol_original.random();
                        sol_original.normalize( mass );
                        
                        const auto& function_rhs  = experiment_rhs;
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        FloatVector rhs_original = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        // rhs_original.zero();
                        
                        // const Float desired_precision = sqrt( machine_epsilon );

                        if( do_cgmpp )
                        {
                            LOG << "CGM C++" << endl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            ConjugateGradientMethod Solver( stiffness );
                            Solver.print_modulo        = 0;
                            Solver.threshold           = desired_precision;
                            Solver.max_iteration_count = 1 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_crmpp_expl )
                        {
                            LOG << "CRM C++" << endl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            ConjugateResidualMethod Solver( stiffness );
                            // Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            // Solver.print_modulo        = 1;
                            Solver.print_modulo        = 0;
                            Solver.threshold           = desired_precision;
                            Solver.max_iteration_count = 1 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve_explicit( sol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_crmpp_robt )
                        {
                            LOG << "CRM C++" << endl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            ConjugateResidualMethod Solver( stiffness );
                            // Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            // Solver.print_modulo        = 1;
                            Solver.print_modulo        = 0;
                            Solver.threshold           = desired_precision;
                            Solver.max_iteration_count = 1 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve_robust( sol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_crmpp_fast )
                        {
                            LOG << "CRM C++" << endl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            ConjugateResidualMethod Solver( stiffness );
                            // Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            // Solver.print_modulo        = 1;
                            Solver.print_modulo        = 0;
                            Solver.threshold           = desired_precision;
                            Solver.max_iteration_count = 1 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve_fast( sol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_minres )
                        {
                            LOG << "MINRES C++" << endl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            MinimumResidualMethod Solver( stiffness );
                            // Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            // Solver.print_modulo        = 1;
                            Solver.print_modulo        = 0;
                            Solver.threshold           = desired_precision;
                            Solver.max_iteration_count = 1 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_herzog )
                        {
                            LOG << "HERZOG SOODHALTER C++" << endl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            HerzogSoodhalterMethod Solver( stiffness );
                            // Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            // Solver.print_modulo        = 1;
                            Solver.print_modulo        = 0;
                            Solver.threshold           = desired_precision;
                            Solver.max_iteration_count = 1 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }
                        
                        if( do_cgm_csr )
                        {
                            LOG << "CGM - CSR Classic" << endl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            FloatVector residual( rhs );
                            auto max_iteration_count = sol.getdimension();
                            timestamp start = gettimestamp();
                            auto recent_iteration_count = 
                            ConjugateGradientSolverCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count ) / max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_crm_csr )
                        {
                            LOG << "CRM - CSR Classic" << endl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            FloatVector residual( rhs );
                            auto max_iteration_count = sol.getdimension();
                            timestamp start = gettimestamp();
                            auto recent_iteration_count = 
                            ConjugateResidualSolverCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count ) / max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_crm_csrtextbook )
                        {
                            LOG << "CRM - CSR Textbook" << endl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            FloatVector residual( rhs );
                            auto max_iteration_count = sol.getdimension();
                            timestamp start = gettimestamp();
                            auto recent_iteration_count = 
                            ConjugateResidualSolverCSR_textbook( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count ) / max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_minres_csr )
                        {
                            LOG << "MINRES CSR" << endl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            FloatVector residual( rhs );
                            auto max_iteration_count = sol.getdimension();
                            timestamp start = gettimestamp();
                            auto recent_iteration_count = 
                            MINRESCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count ) / max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }


                        if( do_whatever_csr )
                        {
                            LOG << "WHATEVER CSR" << endl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            FloatVector residual( rhs );
                            auto max_iteration_count = sol.getdimension();
                            timestamp start = gettimestamp();
                            auto recent_iteration_count = 
                            WHATEVER( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count ) / max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_cgm_diagonal_csr )
                        {
                            LOG << "CGM - CSR Classic with diagonal preconditioning" << endl;
                            
                            auto precon = InverseDiagonalPreconditioner( stiffness );

                            
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            FloatVector residual( rhs );
                            auto max_iteration_count = sol.getdimension();
                            timestamp start = gettimestamp();
                            auto recent_iteration_count = 
                            ConjugateGradientSolverCSR_DiagonalPreconditioner( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                precon.getdiagonal().raw()
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count ) / max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_cgm_ssor_csr )
                        {
                            LOG << "CGM - CSR Classic with SSOR" << endl;
                            
                            auto diagonal = stiffness.diagonal();

                            
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            FloatVector residual( rhs );
                            auto max_iteration_count = sol.getdimension();
                            timestamp start = gettimestamp();
                            auto recent_iteration_count = 
                            ConjugateGradientSolverCSR_SSOR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                diagonal.raw(),
                                1.0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count ) / max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        
                        // contable_sol << nl;
                        contable_res << nl;
                        contable_num << nl;
                        contable_sec << nl;

                        // contable_sol.lg( false );
                        contable_res.lg( false );
                        contable_num.lg( false );
                        contable_sec.lg( false );

                    }
                    
                }

                if( l != max_l ){ 
                    LOG << "Refinement..." << endl;
                    M.uniformrefinement();
                }

                
            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
