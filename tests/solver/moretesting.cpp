

/**/

#include <iostream>
#include <ctime>

#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/simpleoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../solver/iterativesolver.hpp"
// #include "../../solver/crm.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "More testing" );

int main()
{
	LOG << "Unit Test: " << TestName << endl;
	
	{
		LOG << "First Something Simple" << endl;
		
		ScalingOperator S( 10, 3.141 );
		ConjugateResidualMethod CRM(S);

		FloatVector rhs(10), x(10);
		x.random(); rhs.zero();
		
		LOG << x << endl;
		CRM.solve( x, rhs );
	}
	

	{
		LOG << "Now something more complicated." << endl;
		
		int dimension = 100;
		
		FloatVector x( dimension );
		for( int p = 0; p < dimension; p++ )
			x.setentry( p, 3. + p * 5. );
		
		SparseMatrix M( dimension, dimension );
        M.reserve( 3 * dimension );
        
		for( int i = 0; i < dimension; i++ ){
			if( i-1 >= 0 ) 
				M.addentry( i, i-1, 1.25 );
			if( i+1 < dimension ) 
				M.addentry( i, i+1, 1.25 );
			M.addentry( i, i, 2.501 );
		}
		
		M.sortentries();
		
		FloatVector b = M * x;
		
		FloatVector y( dimension );
		y.random();
		
		// LOG << M << endl;
		
		ConjugateResidualMethod CRM(M);
	        CRM.max_iteration_count = 100;
		CRM.threshold = 1e-16;
                
                timestamp start = gettimestamp();
                CRM.solve(y,b);
                timestamp end   = gettimestamp();
                LOG << "time elapsed: " << timestamp2measurement( end - start ) << endl;
                
	}
	
	LOG << "Finished Unit Test: " << TestName << endl;

	return 0;
}
