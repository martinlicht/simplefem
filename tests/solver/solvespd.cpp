

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Solve SPD system: CGM, CRM, MINRES, HerzoogSoodhalter" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
        
        LOG << std::setprecision(5);

        if(true){

            ConvergenceTable contable;
            
            bool do_cgmpp      = false;
            bool do_crmpp_expl = false;
            bool do_crmpp_robt = false;
            bool do_crmpp_fast = false;
            bool do_minres     = false;
            bool do_herzog     = false;
            
            do_cgmpp      = true;
            do_crmpp_expl = true;
            do_crmpp_robt = true;
            do_crmpp_fast = true;
            do_minres     = true;
            do_herzog     = true;
            
            
            
            
            contable << "Index";
            if( do_cgmpp      ) contable << "CGM++"         << "time";
            if( do_crmpp_expl ) contable << "CRM++(expl)"   << "time";
            if( do_crmpp_robt ) contable << "CRM++(robt)"   << "time";
            if( do_crmpp_fast ) contable << "CRM++(fast)"   << "time";
            if( do_minres     ) contable << "MINRES"        << "time";
            if( do_herzog     ) contable << "HERZOG"        << "time";
                     ;
            
            const std::vector<int> Ns = { 16, 32 };

            for( const int N : Ns ){
                
                LOG << "Level: " << N << std::endl;
                
                {
                    
                    LOG << "...assemble matrix" << endl;

                    std::vector< SparseMatrix::MatrixEntry > entries;

                    for( int e = 0; e < N*N; e++ )
                    {
                        int x = e / N;
                        int y = e % N;
                        assert( e == x * N + y );

                        entries.push_back({ x * N + y, x * N + y, 4. / N });

                        if( x != 0   ) entries.push_back({ x * N + y, (x-1) * N + y,   -1. / N });
                        if( x != N-1 ) entries.push_back({ x * N + y, (x+1) * N + y,   -1. / N });
                        if( y != 0   ) entries.push_back({ x * N + y, (x  ) * N + y-1, -1. / N });
                        if( y != N-1 ) entries.push_back({ x * N + y, (x  ) * N + y+1, -1. / N });
                        
                    }

                    auto system = SparseMatrix( N*N, N*N, entries );


                    LOG << "...create solutions and right-hand sides" << endl;

                    const int T = 10;

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
                        
                        contable << static_cast<Float>(N);

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
                            contable << stat_sol << stat_num;
                        }

                        if( do_crmpp_expl )
                        {
                            LOG << "CRM C++ (explicit)" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve_explicit( mysol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable << stat_sol << stat_num;
                        }

                        if( do_crmpp_robt )
                        {
                            LOG << "CRM C++ (robust)" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve_robust( mysol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable << stat_sol << stat_num;
                        }

                        if( do_crmpp_fast )
                        {
                            LOG << "CRM C++ (fast)" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve_fast( mysol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable << stat_sol << stat_num;
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
                            contable << stat_sol << stat_num;
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
                            contable << stat_sol << stat_num;
                        }
                        
                        
                        contable << nl;
                        
                    }
                    
                    contable.lg( false );

                    }

                
                
            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
