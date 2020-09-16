

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"


using namespace std;

int main()
{
        cout << std::unitbuf;
        cout << "Unit Test for Vector class" << endl;
        
        if(true)
        {
        
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
            
            cout << "Add the middle slice:" << endl;
            a.addslice( 1, t.getslice(1,3), 1000. );
            cout << a << endl;
            
            FloatVector e(0);
            cout << "Should be the zero-dimensional vector:" << endl;
            cout << e << endl;
            
        }
        
        if(true)
        {
        
            FloatVector a { 1, 3, 0 };
            FloatVector b { -1, -3, 0 };
            FloatVector c { 2, 4, 1 };
            FloatVector d { -5, -4, -3 };
            FloatVector e { 0,0,0 };
            
            cout << "positive:     (no ) " << a.ispositive() << endl;
            cout << "negative:     (no ) " << a.isnegative() << endl;
            cout << "non-negative: (yes) " << a.isnonnegative() << endl;
            cout << "non-positive: (no ) " << a.isnonpositive() << endl;
            cout << "zero:         (no ) " << a.iszero() << endl;
            
            cout << "positive:     (no ) " << b.ispositive() << endl;
            cout << "negative:     (no ) " << b.isnegative() << endl;
            cout << "non-negative: (no ) " << b.isnonnegative() << endl;
            cout << "non-positive: (yes) " << b.isnonpositive() << endl;
            cout << "zero:         (no ) " << b.iszero() << endl;
            
            cout << "positive:     (yes) " << c.ispositive() << endl;
            cout << "negative:     (no ) " << c.isnegative() << endl;
            cout << "non-positive: (no ) " << c.isnonpositive() << endl;
            cout << "non-negative: (yes) " << c.isnonnegative() << endl;
            cout << "zero:         (no ) " << c.iszero() << endl;
            
            cout << "positive:     (no ) " << d.ispositive()    << endl;
            cout << "negative:     (yes) " << d.isnegative() << endl;
            cout << "non-positive: (yes) " << d.isnonpositive() << endl;
            cout << "non-negative: (no ) " << d.isnonnegative() << endl;
            cout << "zero:         (no ) " << d.iszero() << endl;
            
            cout << "positive:     (no ) " << e.ispositive()    << endl;
            cout << "negative:     (no ) " << e.isnegative()    << endl;
            cout << "non-positive: (yes) " << e.isnonpositive() << endl;
            cout << "non-negative: (yes) " << e.isnonnegative() << endl;
            cout << "zero:         (yes) " << e.iszero()        << endl;
            
            cout << "norm: (3.162)" << a.norm() << space << a.normalize().norm() << endl;
            cout << "norm: (3.162)" << b.norm() << space << a.normalize().norm() << endl;
            cout << "norm: (4.582)" << c.norm() << space << a.normalize().norm() << endl;
            cout << "norm: (7.071)" << d.norm() << space << a.normalize().norm() << endl;
            cout << "norm: (0.000)" << e.norm() << endl;
            
        }
        
        cout << "Finished Unit Test" << endl;

        return 0;
}
