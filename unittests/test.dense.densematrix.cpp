

/**/

#include <iostream>
#include "../basic.hpp"
#include "../dense/densematrix.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Dense Matrix class" << endl;
	
	DenseMatrix A( 3, 4 );
	A.randommatrix();
	
	cout << A << endl;
	cout << A.transpose() << endl;
	cout << 3 * A << endl;
	
	DenseMatrix B( 3, 4 );
	B.randommatrix();
	
	cout << B << endl;
	cout << A + B << endl;
	
	DenseMatrix I3(3,3);
	I3.unitmatrix();
	DenseMatrix I4(4,4);
	I4.unitmatrix();
	cout << I3 << I4 << endl;
	cout << I3 * A << endl;
	cout << A * I4 << endl;
	auto S5 = 5. * I3;
	cout << S5 << endl;
	cout << S5 * A << endl;
	
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
