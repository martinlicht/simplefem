

/**/

#include <iostream>
#include "../basic.hpp"
#include "floatvector.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Vector class" << endl;
	
	FloatVector a(5);
	
	a.check();
	
	a.zero();
	cout << "Should be zero vector:" << endl;
	cout << a << endl;
	
	for( int i = 0; i < 5; i++ )
		a.setentry( i, i+1 );
	cout << "Should be ascending numbers:" << endl;
	cout << a << endl;
	
	FloatVector b(a);
	cout << "Should be the same again:" << endl;
	cout << b << endl;
	
	cout << "Should be multiples of PI:" << endl;
	cout << 3.141 * a << endl;
	
	cout << "Should be negative of original vector:" << endl;
	cout << -a << endl;
	
	FloatVector t(5);
	
	for( int i = 0; i < 5; i++ )
		t.setentry( i, 3. * i+1 );
	cout << "Should be other ascending numbers:" << endl;
	cout << t << endl;
	
	cout << "Next the sum:" << endl;
	cout << a+t << endl;
	cout << "Then the difference:" << endl;
	cout << a-t << endl;
	
	cout << "Now the scalar product with itself:" << endl;
	cout << a*a << endl;
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
