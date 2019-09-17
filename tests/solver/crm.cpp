

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/scalingoperator.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../solver/crm.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Conjugate Residual Method" << endl;
	
	{
		
		cout << "First Something Simple" << endl;
		
		ScalingOperator S( 10, 3.141 );
		ConjugateResidualMethod CRM(S);

		FloatVector rhs(10), x(10);
		x.random(); rhs.zero();
		
		cout << x << endl;
		cout << CRM << endl;
		CRM.solve( x, rhs );
		// cout << x << endl;
	
	}
	

	if(true){
		
		cout << "Now something more complicated." << endl;
		
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
			M.addentry( i, i, 2.51 );
		}
		M.sortentries();
		
		cout << "Compute stuff." << endl;
		
		FloatVector b = M * x;
		
		FloatVector y( dimension );
		y.random();
		
		// cout << M << endl;
		
		ConjugateResidualMethod CRM(M);
		CRM.max_iteration_count = 50;
		CRM.tolerance = 1e-20;
		
		timestamp start, end;
		start = gettimestamp();
		CRM.solve(y,b);
		end = gettimestamp();
		cout << end - start << endl;

		
	}
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
