

/**/

#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/simpleoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../solver/iterativesolver.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Tests: some basic tests of solvers" << nl;

    {
        
        LOG << "First Something Simple..." << nl;
        
        ScalingOperator A( 10, Constants::pi );
        ScalingOperator M( 10, Constants::euler );
        
        FloatVector rhs(10), x(10);

        x.random(); rhs.zero();
        LOG << x << nl;
        ConjugateResidualMethod CRM(A);
        CRM.solve( x, rhs );
        // LOG << x << nl;
        
        x.random(); rhs.zero();
        LOG << x << nl;
        MinimumResidualMethod MINRES(A);
        MINRES.solve( x, rhs );
        // LOG << x << nl;
        
        x.random(); rhs.zero();
        LOG << x << nl;
        PreconditionedConjugateResidualMethod PCRM(A,M);
        PCRM.solve( x, rhs );
        // LOG << x << nl;
        
    }
    
    {
        
        LOG << "Now something more complicated:\n Tridiagonal with very weak diagonal dominance." << nl;
        
        int dimension = 100;
        
        FloatVector x( dimension );
        for( int p = 0; p < dimension; p++ )
            x.setentry( p, 3. + p * 5. );
        
        SparseMatrix A( dimension, dimension );
        A.reserve( 3 * dimension );
        
        for( int i = 0; i < dimension; i++ ){
            if( i-1 >= 0 ) 
                A.appendentry( i, i-1, 1.25 );
            if( i+1 < dimension ) 
                A.appendentry( i, i+1, 1.25 );
            A.appendentry( i, i, 2.51 );
        }
        // A.sortentries();
        
        auto M = IdentityOperator( dimension );
        
        // SparseMatrix M( dimension, dimension );
        // M.reserve( dimension );
        // for( int i = 0; i < dimension; i++ ){
        //     M.appendentry( i, i, 1./2.51 );
        // }
        // // M.sortentries();
        
        LOG << "Compute stuff." << nl;
        
        FloatVector b = A * x;

        FloatVector y_original( dimension ); y_original.random();
        
        {
            FloatVector y = y_original;

            ConjugateGradientMethod CGM(A);
            CGM.verbosity = IterativeSolver::VerbosityLevel::startandfinish;
            CGM.print_modulo = 0;
            CGM.max_iteration_count = 2 * dimension;
            CGM.tolerance = desired_precision;
            
            timestamp start, end;
            start = timestampnow();
            CGM.solve(y,b);
            end = timestampnow();
            LOG << timestamp2measurement( end - start ) << nl;
        }
        
        {
            FloatVector y = y_original;

            PreconditionedConjugateGradientMethod PCGM(A,M);
            PCGM.verbosity = IterativeSolver::VerbosityLevel::startandfinish;
            PCGM.print_modulo = 0;
            PCGM.max_iteration_count = 2 * dimension;
            PCGM.tolerance = desired_precision;
            
            timestamp start, end;
            start = timestampnow();
            PCGM.solve(y,b);
            end = timestampnow();
            LOG << timestamp2measurement( end - start ) << nl;
        }
        
        {
            FloatVector y = y_original;
            
            ConjugateResidualMethod CRM(A);
            CRM.verbosity = IterativeSolver::VerbosityLevel::startandfinish;
            CRM.print_modulo = 0;
            CRM.max_iteration_count = 2 * dimension;
            CRM.tolerance = desired_precision;
            
            timestamp start, end;
            start = timestampnow();
            CRM.solve_explicit(y,b);
            end = timestampnow();
            LOG << timestamp2measurement( end - start ) << nl;
        }
        
        {
            FloatVector y = y_original;
            
            PreconditionedConjugateResidualMethod PCRM(A,M);
            PCRM.verbosity = IterativeSolver::VerbosityLevel::startandfinish;
            PCRM.print_modulo = 0;
            PCRM.max_iteration_count = 2 * dimension;
            PCRM.tolerance = desired_precision;
            
            timestamp start, end;
            start = timestampnow();
            PCRM.solve(y,b);
            end = timestampnow();
            LOG << timestamp2measurement( end - start ) << nl;
        }
        
        {
            FloatVector y = y_original;
            
            MinimumResidualMethod MINRES(A);
            MINRES.verbosity = IterativeSolver::VerbosityLevel::startandfinish;
            MINRES.print_modulo = 0;
            MINRES.max_iteration_count = 2 * dimension;
            MINRES.tolerance = desired_precision;
            
            timestamp start, end;
            start = timestampnow();
            MINRES.solve(y,b);
            end = timestampnow();
            LOG << timestamp2measurement( end - start ) << nl;
        }
        
    }

    {
        
        LOG << "For MINRES: diagonal indefinite matrix ." << nl;
        
        int dimension = 100;
        
        FloatVector x( dimension );
        for( int p = 0; p < dimension; p++ )
            x.setentry( p, 3. + p * 5. );
        
        SparseMatrix A( dimension, dimension );
        A.reserve( dimension );
        
        // for( int i = 0; i < dimension; i++ ){
        //     A.appendentry( i, i, sign_power(i) * 2.51 );
        // }
        
        for( int i = 0; i < dimension/2; i++ ){
            A.appendentry(               i,               i,  3+i );
            A.appendentry( dimension/2 + i, dimension/2 + i, -3-i );
        }
        
        LOG << "Compute stuff." << nl;
        
        FloatVector b = A * x;
        
        {
            FloatVector y( dimension ); srand(0); y.random();
            
            MinimumResidualMethod MINRES(A);
            MINRES.max_iteration_count = 2 * dimension;
            MINRES.tolerance = desired_precision;
            
            timestamp start, end;
            start = timestampnow();
            MINRES.solve(y,b);
            end = timestampnow();
            LOG << timestamp2measurement( end - start ) << nl;
        }
        
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
