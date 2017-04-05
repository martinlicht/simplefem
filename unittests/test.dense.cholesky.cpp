

/**/

#include <iostream>
#include "../basic.hpp"
#include "../dense/dense.functions.hpp"
#include "../dense/qr.factorization.hpp"
#include "../dense/lu.factorization.hpp"
#include "../dense/cholesky.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Cholesky decomposition" << endl;
	
	std::cout.precision(5);
	std::cout.setf( std::ios::fixed, std:: ios::floatfield );
	std::cout << std::showpos;
	
	DenseMatrix A(3,3);
	
        A(0,0) =   4; A(0,1) =  12; A(0,2) = -16; 
        A(1,0) =  12; A(1,1) =  37; A(1,2) = -43; 
        A(2,0) = -16; A(2,1) = -43; A(2,2) =  98; 
        
        DenseMatrix L = CholeskyDecomposition( A );
        
        cout << "Original matrix:" << A << endl;
        
        cout << "Factor matrix:" << L << endl;
        
        cout << "Product matrix:" << L * Transpose(L) << endl;
        
        cout << "Finished Unit Test" << endl;

	return 0;
}
