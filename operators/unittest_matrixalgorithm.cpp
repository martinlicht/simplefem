

/**/

#include <iostream>
#include "../basic.hpp"
#include "matrixalgorithm.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Matrix Algorithms class" << endl;
	
	// {
		// DenseMatrix A(3,3);
		// for( int t = 0; t < 3; t++ )
		// for( int s = 0; s < 3; s++ )
			// A(s,t) = 2 + s*s + t*t*t;
		
		// DenseMatrix Ainv = Inverse( A );
		
		// cout << A << Ainv << A * Ainv << endl;
	// }
	
	{
		const int D = 10;
		DenseMatrix A(D,D);
		for( int t = 0; t < D; t++ )
		for( int s = 0; s < D; s++ )
			A(s,t) = 2 + s*s + t*t*t;
		
		DenseMatrix Q(D,D), R(D,D);
		
		PolarDecomposition( A, Q, R );
		
		cout << A << Q << R << Q * R << endl;
	}
	
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
