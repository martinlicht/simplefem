

/**/

#include <ostream>
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/simpleoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../solver/iterativesolver.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/pcrm.hpp"


using namespace std;

int main()
{
    LOG << "Unit Tests: some basic tests of Solvers" << endl;

    {
        
        LOG << "First Something Simple..." << endl;
        
        ScalingOperator A( 10, Constants::pi );
        ScalingOperator M( 10, Constants::euler );
        
        FloatVector rhs(10), x(10);

        x.random(); rhs.zero();
        LOG << x << endl;
        ConjugateResidualMethod CRM(A);
        CRM.solve( x, rhs );
        // LOG << x << endl;
        
        x.random(); rhs.zero();
        LOG << x << endl;
        MinimumResidualMethod MINRES(A);
        MINRES.solve( x, rhs );
        // LOG << x << endl;
        
        x.random(); rhs.zero();
        LOG << x << endl;
        PreconditionedConjugateResidualMethod PCRM(A,M);
        PCRM.solve( x, rhs );
        // LOG << x << endl;
        
    }
    
    {
        
        LOG << "Now something more complicated:\n Tridiagonal with very weak diagonal dominance." << endl;
        
        int dimension = 100;
        
        FloatVector x( dimension );
        for( int p = 0; p < dimension; p++ )
            x.setentry( p, 3. + p * 5. );
        
        SparseMatrix A( dimension, dimension );
        A.reserve( 3 * dimension );
        
        for( int i = 0; i < dimension; i++ ){
            if( i-1 >= 0 ) 
                A.addentry( i, i-1, 1.25 );
            if( i+1 < dimension ) 
                A.addentry( i, i+1, 1.25 );
            A.addentry( i, i, 2.51 );
        }
        // A.sortentries();
        
        SparseMatrix M( dimension, dimension );
        M.reserve( dimension );
        
        for( int i = 0; i < dimension; i++ ){
            M.addentry( i, i, 1./2.51 );
        }
        // A.sortentries();
        
        LOG << "Compute stuff." << endl;
        
        FloatVector b = A * x;
        
        {
            FloatVector y( dimension );
            srand(0);
            y.random();
            
            ConjugateResidualMethod CRM(A);
            CRM.max_iteration_count = dimension;
            CRM.threshold = 1e-20;
            
            timestamp start, end;
            start = gettimestamp();
            CRM.solve(y,b);
            end = gettimestamp();
            LOG << timestamp2measurement( end - start ) << endl;
        }
        
        {
            FloatVector y( dimension );
            srand(0);
            y.random();
            
            PreconditionedConjugateResidualMethod PCRM(A,M);
            PCRM.max_iteration_count = dimension;
            PCRM.threshold = 1e-20;
            
            timestamp start, end;
            start = gettimestamp();
            PCRM.solve(y,b);
            end = gettimestamp();
            LOG << timestamp2measurement( end - start ) << endl;
        }
        
        {
            FloatVector y( dimension );
            srand(0);
            y.random();
            
            MinimumResidualMethod MINRES(A);
            MINRES.max_iteration_count = dimension;
            MINRES.threshold = 1e-20;
            
            timestamp start, end;
            start = gettimestamp();
            MINRES.solve(y,b);
            end = gettimestamp();
            LOG << timestamp2measurement( end - start ) << endl;
        }
        
    }

    {
        
        LOG << "For MINRES: diagonal indefinite matrix ." << endl;
        
        int dimension = 100;
        
        FloatVector x( dimension );
        for( int p = 0; p < dimension; p++ )
            x.setentry( p, 3. + p * 5. );
        
        SparseMatrix A( dimension, dimension );
        A.reserve( dimension );
        
        // for( int i = 0; i < dimension; i++ ){
        //     A.addentry( i, i, signpower(i) * 2.51 );
        // }
        
        for( int i = 0; i < dimension/2; i++ ){
            A.addentry(               i,               i,  3+i );
            A.addentry( dimension/2 + i, dimension/2 + i, -3-i );
        }
        
        LOG << "Compute stuff." << endl;
        
        FloatVector b = A * x;
        
        {
            FloatVector y( dimension );
            srand(0);
            y.random();
            
            MinimumResidualMethod MINRES(A);
            MINRES.max_iteration_count = dimension;
            MINRES.threshold = 1;
            
            timestamp start, end;
            start = gettimestamp();
            MINRES.solve(y,b);
            end = gettimestamp();
            LOG << timestamp2measurement( end - start ) << endl;
        }
        
    }
    
    LOG << "Finished Unit Test" << endl;

    return 0;
}
