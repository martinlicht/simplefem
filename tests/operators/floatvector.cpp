

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Float vector class" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
        
        if(true)
        {
        
            FloatVector a(5);
            
            a.check();
            
            LOG << "Should be zero vector:" << endl;
            a.zero();
            LOG << a << endl;
            
            for( int i = 0; i < 5; i++ )
                    a.setentry( i, i+1 );
            LOG << "Should be ascending numbers:" << endl;
            LOG << a << endl;
            
            LOG << "Should be the middle entries:" << endl;
            LOG << a.getslice(1,3) << endl;
            
            FloatVector b(a);
            LOG << "Should be the same again:" << endl;
            LOG << b << endl;
            
            LOG << "Should be multiples of PI:" << endl;
            LOG << 3.141 * a << endl;
            
            LOG << "Should be negative of original vector:" << endl;
            LOG << -a << endl;
            
            FloatVector t(5);
            
            for( int i = 0; i < 5; i++ )
                    t.setentry( i, 3. * i+1 );
            LOG << "Should be other ascending numbers: 3*( i + 1) " << endl;
            LOG << t << endl;
            
            LOG << "Next the sum of two vectors:" << endl;
            LOG << a+t << endl;
            LOG << "Then their difference:" << endl;
            LOG << a-t << endl;
            
            LOG << "Now the scalar product with itself:" << endl;
            LOG << a*a << endl;
            
            LOG << "Copy the middle slice:" << endl;
            a.setslice( 1, t.getslice(1,3) );
            LOG << a << endl;
            
            LOG << "Add the middle slice:" << endl;
            a.addslice( 1, t.getslice(1,3), 1000. );
            LOG << a << endl;
            
            FloatVector e(0);
            LOG << "Should be the zero-dimensional vector:" << endl;
            LOG << e << endl;
            
        }
        
        if(true)
        {
        
            FloatVector a { 1, 3, 0 };
            FloatVector b { -1, -3, 0 };
            FloatVector c { 2, 4, 1 };
            FloatVector d { -5, -4, -3 };
            FloatVector e { 0,0,0 };
            
            LOG << "positive:     (no ) " << a.ispositive() << endl;
            LOG << "negative:     (no ) " << a.isnegative() << endl;
            LOG << "non-negative: (yes) " << a.isnonnegative() << endl;
            LOG << "non-positive: (no ) " << a.isnonpositive() << endl;
            LOG << "zero:         (no ) " << a.iszero() << endl;
            
            LOG << "positive:     (no ) " << b.ispositive() << endl;
            LOG << "negative:     (no ) " << b.isnegative() << endl;
            LOG << "non-negative: (no ) " << b.isnonnegative() << endl;
            LOG << "non-positive: (yes) " << b.isnonpositive() << endl;
            LOG << "zero:         (no ) " << b.iszero() << endl;
            
            LOG << "positive:     (yes) " << c.ispositive() << endl;
            LOG << "negative:     (no ) " << c.isnegative() << endl;
            LOG << "non-positive: (no ) " << c.isnonpositive() << endl;
            LOG << "non-negative: (yes) " << c.isnonnegative() << endl;
            LOG << "zero:         (no ) " << c.iszero() << endl;
            
            LOG << "positive:     (no ) " << d.ispositive()   << endl;
            LOG << "negative:     (yes) " << d.isnegative() << endl;
            LOG << "non-positive: (yes) " << d.isnonpositive() << endl;
            LOG << "non-negative: (no ) " << d.isnonnegative() << endl;
            LOG << "zero:         (no ) " << d.iszero() << endl;
            
            LOG << "positive:     (no ) " << e.ispositive()   << endl;
            LOG << "negative:     (no ) " << e.isnegative()   << endl;
            LOG << "non-positive: (yes) " << e.isnonpositive() << endl;
            LOG << "non-negative: (yes) " << e.isnonnegative() << endl;
            LOG << "zero:         (yes) " << e.iszero()       << endl;
            
            LOG << "norm: (3.162)" << a.norm() << space << a.normalize().norm() << endl;
            LOG << "norm: (3.162)" << b.norm() << space << a.normalize().norm() << endl;
            LOG << "norm: (4.582)" << c.norm() << space << a.normalize().norm() << endl;
            LOG << "norm: (7.071)" << d.norm() << space << a.normalize().norm() << endl;
            LOG << "norm: (0.000)" << e.norm() << endl;
            
        }
        
        std::clog << "Finished Unit Test: " << TestName << endl;

        return 0;
}
