

/**/

#include <iostream>
#include "../basic.hpp"
#include "../operators/floatvector.hpp"


using namespace std;

int main()
{
        cout << std::unitbuf;
        cout << "Unit Test for Vector class" << endl;
        
        FloatVector a(5);
        
        a.check();
        
        cout << "Should be zero vector:" << endl;
        a.zero();
        cout << a << endl;
        
        cout << "Print plain:" << endl;
        a.printplain( cout );
        
        for( int i = 0; i < 5; i++ )
                a.setentry( i, i+1 );
        cout << "Should be ascending numbers:" << endl;
        cout << a << endl;
        
        cout << "Should be the middle entries:" << endl;
        cout << a.getslice(1,3) << endl;
        
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
        cout << "Should be other ascending numbers: 3*( i + 1) " << endl;
        cout << t << endl;
        
        cout << "Next the sum of two vectors:" << endl;
        cout << a+t << endl;
        cout << "Then their difference:" << endl;
        cout << a-t << endl;
        
        cout << "Now the scalar product with itself:" << endl;
        cout << a*a << endl;
        
        cout << "Copy the middle slice:" << endl;
        a.setslice( 1, t.getslice(1,3) );
        cout << a << endl;
        
        FloatVector e(0);
        cout << "Should be the zero-dimensional vector:" << endl;
        cout << e << endl;
        
        cout << "Finished Unit Test" << endl;

        return 0;
}
