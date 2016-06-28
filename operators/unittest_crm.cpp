

/**/

#include <iostream>
#include "../basic.hpp"
#include "floatvector.hpp"
#include "scalingoperator.hpp"
#include "crm.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Conjugate Residual Method" << endl;
	
	ScalingOperator S( 10, 3.141 );
	ConjugateResidualMethod CRM(S);

	FloatVector rhs(10), x(10);
	x.random(); rhs.zero();
	
	cout << x << endl;
	
	CRM.solve( x, rhs );
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
