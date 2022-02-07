

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../solver/chebyshev.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Solve SPD system: CSR Solvers" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
        
        LOG << std::setprecision(5);

        if(true){

            ConvergenceTable contable_sol("L2 Error");
            ConvergenceTable contable_res("L2 Residual");
            ConvergenceTable contable_num("Iteration percentage");
            
            contable_sol.print_transpose_instead_of_standard = true;
            contable_res.print_transpose_instead_of_standard = true;
            contable_num.print_transpose_instead_of_standard = true;
            
            bool do_cgmpp      = false;
            bool do_crmpp_expl = false;
            bool do_crmpp_robt = false;
            bool do_crmpp_fast = false;
            bool do_minres     = false;
            bool do_herzog     = false;
            //
            bool do_cgm_csr                = false;
            bool do_crm_csr                = false;
            bool do_crm_csrtextbook        = false;
            bool do_minres_csr             = false;
            bool do_whatever_csr           = false;
            bool do_cgm_diagonal_csr       = false;
            bool do_cgm_ssor_csr           = false;
            bool do_chebyshev_diagonal_csr = false;

            do_cgmpp      = true;
            do_crmpp_expl = true;
            do_crmpp_robt = true;
            do_crmpp_fast = true;
            do_minres     = true;
            do_herzog     = true;
            //
            do_cgm_csr                = true;
            do_crm_csr                = true;
            do_crm_csrtextbook        = true;
            do_minres_csr             = true;
            do_whatever_csr           = true;
            do_cgm_diagonal_csr       = true;
            do_cgm_ssor_csr           = true;
            do_chebyshev_diagonal_csr = true;
            
            
            
            
            contable_sol << "Index";
            if( do_cgmpp      ) contable_sol << "CGM++"      ;
            if( do_crmpp_expl ) contable_sol << "CRM++(expl)";
            if( do_crmpp_robt ) contable_sol << "CRM++(robt)";
            if( do_crmpp_fast ) contable_sol << "CRM++(fast)";
            if( do_minres     ) contable_sol << "MINRES"     ;
            if( do_herzog     ) contable_sol << "HERZOG"     ;
            //
            if( do_cgm_csr )                contable_sol << "CGMcsr"       ;
            if( do_crm_csr )                contable_sol << "CRMcsr"       ;
            if( do_crm_csrtextbook )        contable_sol << "CRMcsr_tb"    ;
            if( do_minres_csr )             contable_sol << "MINREScsr"    ;
            if( do_whatever_csr )           contable_sol << "WHATEVER"     ;
            if( do_cgm_diagonal_csr )       contable_sol << "CGMcsr_diag"  ;
            if( do_cgm_ssor_csr )           contable_sol << "CGMcsr_ssor"  ;
            if( do_chebyshev_diagonal_csr ) contable_sol << "Chebyshev_csr";
            
            contable_res << "Index";
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

            contable_num << "Index";
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
            

            const std::vector<int> Ns = { 4, 8, 16, 32 };

            for( const int N : Ns ){
                
                LOG << "Level: " << N << std::endl;
                
                {
                    
                    LOG << "...assemble matrix" << endl;

                    std::vector< SparseMatrix::MatrixEntry > entries;

                    const Float h = 1. / (N+1);
                    const Float h2 = h * h;

                    for( int e = 0; e < N*N; e++ )
                    {
                        int x = e / N;
                        int y = e % N;
                        assert( e == x * N + y );

                        entries.push_back({ x * N + y, x * N + y, 4. / h2 });

                        if( x != 0   ) entries.push_back({ x * N + y, (x-1) * N + y,   -1. / h2 });
                        if( x != N-1 ) entries.push_back({ x * N + y, (x+1) * N + y,   -1. / h2 });
                        if( y != 0   ) entries.push_back({ x * N + y, (x  ) * N + y-1, -1. / h2 });
                        if( y != N-1 ) entries.push_back({ x * N + y, (x  ) * N + y+1, -1. / h2 });
                        
                    }

                    auto system_prelim = SparseMatrix( N*N, N*N, entries );
                    system_prelim.sortentries();
                    auto system = MatrixCSR( system_prelim );


                    LOG << "...create solutions and right-hand sides" << endl;

                    const int T = 5;

                    std::vector<FloatVector> sols;
                    std::vector<FloatVector> rhss;

                    for( int t = 0; t < T; t++ )
                    {
                        FloatVector sol( N * N );
                        sol.random();
                        sol.normalize();
                        FloatVector rhs = system * sol;
                        sols.push_back( sol );
                        rhss.push_back( rhs );
                    }



                    for( int t = 0; t < T; t++ )
                    {

                        const auto& sol = sols[t];
                        const auto& rhs = rhss[t];
                        
                        contable_sol << static_cast<Float>(N);
                        contable_res << static_cast<Float>(N);
                        contable_num << static_cast<Float>(N);

                        if( do_cgmpp )
                        {
                            LOG << "CGM C++" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateGradientMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( mysol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }

                        if( do_crmpp_expl )
                        {
                            LOG << "CRM C++ (explicit)" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count = 4 * mysol.getdimension();
                            
                            timestamp start = gettimestamp();
                            Solver.solve_explicit( mysol, rhs );
                            timestamp end = gettimestamp();
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }

                        if( do_crmpp_robt )
                        {
                            LOG << "CRM C++ (robust)" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count = 4 * mysol.getdimension();
                            
                            timestamp start = gettimestamp();
                            Solver.solve_robust( mysol, rhs );
                            timestamp end = gettimestamp();
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }

                        if( do_crmpp_fast )
                        {
                            LOG << "CRM C++ (fast)" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count = 4 * mysol.getdimension();
                            
                            timestamp start = gettimestamp();
                            Solver.solve_fast( mysol, rhs );
                            timestamp end = gettimestamp();
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }

                        if( do_minres )
                        {
                            LOG << "MINRES C++" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            MinimumResidualMethod Solver( system );
                            Solver.print_modulo        = 4 * mysol.getdimension();
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( mysol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }

                        if( do_herzog )
                        {
                            LOG << "HERZOG SOODHALTER C++" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            HerzogSoodhalterMethod Solver( system );
                            Solver.print_modulo        = 4 * mysol.getdimension();
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( mysol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }
                        
                        
                        




                        if( do_cgm_csr )
                        {
                            LOG << "CGM - CSR Classic" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = gettimestamp();
                            int recent_iteration_count =
                            ConjugateGradientSolverCSR( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            timestamp end = gettimestamp();

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }

                        if( do_crm_csr )
                        {
                            LOG << "CRM - CSR Classic" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = gettimestamp();
                            int recent_iteration_count =
                            ConjugateResidualSolverCSR( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            timestamp end = gettimestamp();

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }

                        if( do_crm_csrtextbook )
                        {
                            LOG << "CRM - CSR Textbook" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = gettimestamp();
                            int recent_iteration_count =
                            ConjugateResidualSolverCSR_textbook( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            timestamp end = gettimestamp();

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }

                        if( do_minres_csr )
                        {
                            LOG << "MINRES CSR" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = gettimestamp();
                            int recent_iteration_count =
                            MINRESCSR( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            timestamp end = gettimestamp();

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }


                        if( do_whatever_csr )
                        {
                            LOG << "WHATEVER CSR" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = gettimestamp();
                            int recent_iteration_count =
                            WHATEVER( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            timestamp end = gettimestamp();

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }


                        if( do_cgm_diagonal_csr )
                        {
                            LOG << "CGM diagonal preconditioner CSR" << endl;
                        
                            DiagonalOperator invprecon = InverseDiagonalPreconditioner( system_prelim );
//                             invprecon.setentries( 1. );
                            assert( invprecon.getdiagonal().isfinite() );
                            assert( invprecon.getdiagonal().isnonnegative() );
                            
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = gettimestamp();
                            int recent_iteration_count =
                            ConjugateGradientSolverCSR_DiagonalPreconditioner( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                invprecon.getdiagonal().raw()
                            );
                            timestamp end = gettimestamp();

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }
                        
                        
                        if( do_cgm_ssor_csr )
                        {
                            LOG << "CGM SSOR preconditioner CSR" << endl;
                        
                            FloatVector diagonal = system.diagonal();
                            assert( diagonal.isfinite() );
                            assert( diagonal.isnonnegative() );
                            
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = gettimestamp();
                            int recent_iteration_count =
                            ConjugateGradientSolverCSR_SSOR( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                diagonal.raw(),
                                0.9123456789
                            );
                            timestamp end = gettimestamp();

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }
                        
                        
                        if( do_chebyshev_diagonal_csr )
                        {
                            LOG << "CHEBYSHEV CSR" << endl;
                        
                            DiagonalOperator invprecon = InverseDiagonalPreconditioner( system_prelim );
                            assert( invprecon.getdiagonal().isfinite() );
                            assert( invprecon.getdiagonal().ispositive() );
                            
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = gettimestamp();
                            int recent_iteration_count =
                            CheybyshevIteration_DiagonalPreconditioner( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                1,
                                invprecon.getdiagonal().raw(),
                                0.,
                                system.eigenvalueupperbound()
                            );
                            timestamp end = gettimestamp();

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }
                        
                        contable_sol << nl;
                        contable_res << nl;
                        contable_num << nl;
                        
                    }
                    
                    contable_sol.lg( false );
                    contable_res.lg( false );
                    contable_num.lg( false );

                    }

                
                
            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
