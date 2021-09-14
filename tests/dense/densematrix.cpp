

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Dense Matrix class" );

int main()
{
	LOG << "Unit Test: " << TestName << endl;
	
	DenseMatrix A( 3, 4 );
	A.randommatrix();
	
	LOG << A << endl;
	LOG << 3 * A << endl;
	
	DenseMatrix B( 3, 4 );
	B.randommatrix();
	
	LOG << B << endl;
	LOG << A + B << endl;
	
	DenseMatrix I3(3,3);
	I3.unitmatrix();
	DenseMatrix I4(4,4);
	I4.unitmatrix();
	LOG << I3 << I4 << endl;
	LOG << I3 * A << endl;
	LOG << A * I4 << endl;
	auto S5 = 5. * I3;
	LOG << S5 << endl;
	LOG << S5 * A << endl;
	
	
	LOG << "Finished Unit Test: " << TestName << endl;

	return 0;
}
