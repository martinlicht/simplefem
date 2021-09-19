

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

TESTNAME( "Solve SPD system: CGM, CRM, MINRES, HerzoogSoodhalter, CSR Solvers" );

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

                        const auto& rhs = rhss[t];
                        
                        contable << static_cast<Float>(N);

                        if(false)
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

                        if(false)
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

                        if(false)
                        {
                            LOG << "MINRES C++" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            MinimumResidualMethod Solver( system );
                            Solver.print_modulo        = 0;
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
                            Solver.print_modulo        = 0;
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
                            LOG << "CGM - CSR Classic" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
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
                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }

                        if(false)
                        {
                            LOG << "CRM - CSR Classic" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
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
                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }

                        if(false)
                        {
                            LOG << "CRM - CSR Textbook" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
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
                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }

                        if(false)
                        {
                            LOG << "MINRES CSR" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
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
                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }


                        if(false)
                        {
                            LOG << "WHATEVER CSR" << endl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
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
                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }


                        if(false)
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
                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }
                        
                        
                        if(false)
                        {
                            LOG << "CGM SSOR preconditioner CSR" << endl;
                        
                            FloatVector diagonal = system.diagonal();
                            assert( diagonal.isfinite() );
                            assert( diagonal.isnonnegative() );
                            
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
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
                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }
                        
                        
                        if(false)
                        {
                            LOG << "CHEBYSHEV CSR" << endl;
                        
                            DiagonalOperator invprecon = InverseDiagonalPreconditioner( system_prelim );
                            assert( invprecon.getdiagonal().isfinite() );
                            assert( invprecon.getdiagonal().ispositive() );
                            
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            CheybyshevIteration_DiagonalPreconditioner( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                10,
                                invprecon.getdiagonal().raw(),
                                0.,
                                100 * invprecon.getdiagonal().maxnorm()
                            );

                            timestamp end = gettimestamp();
                            auto c1 = static_cast<Float>( end - start );
                            auto c2 = Float( ( system * mysol - rhs ).norm() );
                            contable << c1 << c2;
                        }
                        

                        //if(false)
                        {
                            LOG << "CHEBYSHEV CSR" << endl;
                        
                            DiagonalOperator invprecon = InverseDiagonalPreconditioner( system_prelim );
//                             invprecon.setentries( 1. );
                            assert( invprecon.getdiagonal().isfinite() );
                            assert( invprecon.getdiagonal().isnonnegative() );
                            
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            
                            CHEBY( system, mysol, rhs,
                                IdentityOperator(N*N),
                                10 * mysol.getdimension(), 10-10,
                                0.000000, 8.000001 / h2 );
      
                            // Chebyshev( 
                            //     system, 
                            //     rhs, 
                            //     mysol, 
                            //     residual, 
                            //     10 * mysol.getdimension(), 
                            //     0.00001, 
                            //     18.000001 / h2
                            // );
                            // CheybyshevIteration_DiagonalPreconditioner( 
                            //     mysol.getdimension(), 
                            //     mysol.raw(), 
                            //     rhs.raw(), 
                            //     system.getA(), system.getC(), system.getV(),
                            //     residual.raw(),
                            //     desired_precision,
                            //     10,
                            //     invprecon.getdiagonal().raw(),
                            //     0.,
                            //     18.000001 / h2
                            // );

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
