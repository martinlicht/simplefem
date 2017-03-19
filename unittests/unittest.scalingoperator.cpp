

/**/

#include <iostream>
#include "../basic.hpp"
#include "floatvector.hpp"
#include "scalingoperator.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for ScalingOperator" << endl;
	
	FloatVector a(5);
	for( int i = 0; i < 5; i++ )
		a.setentry( i, i+1 );
	
	ScalingOperator S( 5, 3.14159 );
	
	cout << "We start with this Vector:" << endl;
	cout << a << endl;
	
	cout << "Scaled with PI:" << endl;
	cout << S * a << endl;
	
	cout << "The Scaling is " << S.getscaling() << endl;
	S.setscaling( 2.718 );
	cout << "Now the Snmcaling is " << S.getscaling() << endl;
	cout << "Accordingly:" << endl;
	cout << S * a << endl;
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
