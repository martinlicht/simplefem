

/**/

#include <iostream>
#include "../basic.hpp"
#include "../dense/matrixalgorithm.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Matrix Algorithms class" << endl;
	
	std::cout.precision(5);
	std::cout.setf( std::ios::fixed, std:: ios::floatfield );
	std::cout << std::showpos;
	
	const int D = 8;
	DenseMatrix A(D,D);
	
        A.zeromatrix();
	for( int s = 0; s < D; s++ )
	for( int t = 0; t < D; t++ )
		A(s,t) = 3 * kronecker(s,t) - kronecker(s,t-1) - kronecker(s,t+1);
		
	{
		DenseMatrix Q(D,D), R(D,D);
		
		PolarDecomposition( A, Q, R );
		
		// for( int r = 0; r < D; r++ )
		// for( int c = 0; c < D; c++ )
			// cout << ;
		
		DenseMatrix Rinv = UpperTriangularInverse(R);
		DenseMatrix Qt = Q.transpose();
		
		cout << "Matrix A:" << A;
		cout << "Matrix Q:" << Q;
		cout << "Matrix R:" << R;
		cout << "Matrix Q * R:" << Q * R;
		cout << "Matrix Q^t:" << Qt;
		cout << "Matrix Rinv:" << Rinv;
		cout << "Matrix R * Rinv:" << R * Rinv;
		cout << "Matrix Rinv * R:" << Rinv * R;
		cout << "Matrix Q * Qinv:" << Q * Qt;
		cout << "Matrix Qinv * Q:" << Qt * Q;
		cout << "Matrix Rinv * Q^t * A:" << Rinv * Qt * A;
		cout << "Matrix A * Rinv * Q^t:" << A * Rinv * Qt;

		
	}

	{
		
                cout << endl << "Now for Cholesky" << endl;
                
                DenseMatrix M = A.transpose() * A;
		M = M * M;
		DenseMatrix R = CholeskyDecomposition( M );
		DenseMatrix Rinv = UpperTriangularInverse(R);
		
		cout << "Matrix M:" << M;
		cout << "Matrix R:" << R;
		cout << "Matrix Rt * R:" << R.transpose() * R;
		cout << "Matrix Rinv * Rtinv:" << Rinv * Rinv.transpose();
		cout << "Matrix M * Rinv * Rtinv:" << M * Rinv * Rinv.transpose();
		cout << "Matrix Rinv * Rtinv * M:" << Rinv * Rinv.transpose() * M;
                
	}
	
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
