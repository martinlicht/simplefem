

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

TESTNAME( "Solve SPD system: CGM, CRM, MINRES, HerzoogSoodhalter" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
        
        LOG << std::setprecision(5);

        if(true){

            ConvergenceTable contable;
            
            const std::vector<int> Ns = { 4, 8, 16, 32 };

            for( const int N : Ns ){
                
                LOG << "Level: " << N << std::endl;
                
                {
                    
                    LOG << "...assemble matrix" << endl;

                    std::vector< SparseMatrix::MatrixEntry > entries;

                    const Float h = 1./N;
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

                        const auto& rhs = rhss[t];
                        
                        contable << static_cast<Float>(N);

                        {
                            LOG << "CGM C++" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateGradientMethod Solver( system );
                            Solver.print_modulo        = 0;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( mysol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }

                        {
                            LOG << "CRM C++" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = 0;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( mysol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }

                        {
                            LOG << "MINRES C++" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            MinimumResidualMethod Solver( system );
                            Solver.print_modulo        = 1;
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( mysol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }

                        if(false)
                        {
                            LOG << "HERZOG SOODHALTER C++" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            HerzogSoodhalterMethod Solver( system );
                            Solver.print_modulo        = 1;
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( mysol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }





                        {
                            LOG << "CHEBYSHEV CSR" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            
                            timestamp start = gettimestamp();
                            CHEBY( system, mysol, rhs,
                                IdentityOperator(N*N),
                                10 * mysol.getdimension(), 10-10,
                                0.000001, 8.000001 / h2 );
                            timestamp end = gettimestamp();
      
                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
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
