

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/simpleoperators.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../sparse/sparsematrix.hpp"
// #include "../../solver/minres.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Minimal Residual Method" );

int main()
{
	LOG << "Unit Test: " << TestName << endl;

    {
        
        LOG << "First Something Simple with the MINRES" << endl;
        
        ScalingOperator A( 10, 3.141 );
        MinimumResidualMethod MINRES(A);

        FloatVector rhs(10), x(10);
        x.random(); rhs.zero();
        
        LOG << x << endl;
        MINRES.solve( x, rhs );
        // LOG << x << endl;

    }



    
    if(true){
        
        LOG << "Now something more complicated:\n Tridiagonal with very weak diagonal dominance." << endl;
        
        int dimension = 100;
        
        FloatVector x( dimension );
        for( int p = 0; p < dimension; p++ )
            x.setentry( p, 3. + p * 5. );
        
        SparseMatrix A( dimension, dimension );
        A.reserve( 3 * dimension );
        
        for( int i = 0; i < dimension; i++ ){
            if( i-1 >= 0         ) A.addentry( i, i-1, 1.25 );
            if( i+1 <  dimension ) A.addentry( i, i+1, 1.25 );
            A.addentry( i, i, 2.51 );
        }
        
        LOG << "Compute stuff." << endl;
        
        FloatVector b = A * x;
        
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
    
    
    if(true){
        
        LOG << "For MINRES: diagonal indefinite matrix ." << endl;
        
        int dimension = 100;
        
        FloatVector x( dimension );
        for( int p = 0; p < dimension; p++ )
            x.setentry( p, 3. + p * 5. );
        
        SparseMatrix A( dimension, dimension );
        A.reserve( dimension );
        
//         for( int i = 0; i < dimension; i++ ){
//             A.addentry( i, i, signpower(i) * 2.51 );
//         }
        
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

    LOG << "Finished Unit Test: " << TestName << endl;

    return 0;
}
