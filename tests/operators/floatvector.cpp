
#include <cmath>

#include <vector>

#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test: Vector class" << nl;

    {
        // Test default initialization
        FloatVector v1(5);
        for( int i = 0; i < 5; i++ ) {
            assert( std::isnan( v1.getentry(i) ) );
        }

        // Test initialization with specific value
        FloatVector v2(5, 3.0);
        for( int i = 0; i < 5; i++ ) {
            assert( v2.getentry(i) == 3.0 );
        }

        // Test initialization from std::vector
        std::vector<Float> vals = {1.0, 2.0, 3.0, 4.0, 5.0};
        FloatVector v3(vals);
        for( size_t i = 0; i < vals.size(); i++ ) {
            assert( v3.getentry(i) == vals[i] );
        }

        // Set entries and check values
        FloatVector v4(5, 0.0);
        for( int i = 0; i < 5; i++ ) {
            v4.setentry( i, static_cast<Float>(i) );
            assert( v4.getentry(i) == i );
        }

    }
    

    {

        FloatVector a(5);

        a.check();

        LOG << "Should be zero vector:" << nl;
        a.zero();
        LOG << a << nl;
        for( int i = 0; i < a.getdimension(); i++ ) assert( a[i] == 0. );

        for( int i = 0; i < 5; i++ )
        a.setentry( i, i+1 );
        LOG << "Should be ascending numbers:" << nl;
        LOG << a.data_as_text() << nl;
        for( int i = 1; i < a.getdimension(); i++ ) assert( a[i] > a[i-1] );

        LOG << "Should be the middle entries:" << nl;
        const auto am = a.getslice(1,3);
        LOG << am.data_as_text() << nl;
        assert( am[0] == a[1] and am[1] == a[2] and am[2] == a[3] );

        FloatVector b(a);
        LOG << "Should be the same again:" << nl;
        LOG << b.data_as_text() << nl;
        for( int i = 1; i < a.getdimension(); i++ ) assert( b[i] == a[i] );

        LOG << "Should be multiples of PI:" << nl;
        const auto api = Constants::pi * a;
        LOG << api.data_as_text() << nl;
        for( int i = 1; i < a.getdimension(); i++ ) assert( api[i] == Constants::pi * a[i] );

        LOG << "Should be negative of original vector:" << nl;
        const auto nega = -a;
        LOG << nega.data_as_text() << nl;
        for( int i = 1; i < a.getdimension(); i++ ) assert( nega[i] == - a[i] );

        FloatVector t(5);

        for( int i = 0; i < 5; i++ )
        t.setentry( i, 3. * i+1 );
        LOG << "Should be other ascending numbers: 3*( i + 1) " << nl;
        LOG << t.data_as_text() << nl;
        for( int i = 1; i < t.getdimension(); i++ ) assert( t[i] == 3. * i + 1. );

        LOG << "Next the sum of two vectors:" << nl;
        const auto apt = a+t;
        LOG << apt.data_as_text() << nl;
        for( int i = 1; i < t.getdimension(); i++ ) assert( apt[i] == a[i] + t[i] );

        LOG << "Then their difference:" << nl;
        const auto amt = a-t;
        LOG << amt.data_as_text() << nl;
        for( int i = 1; i < t.getdimension(); i++ ) assert( amt[i] == a[i] - t[i] );

        LOG << "Now the scalar product with itself:" << nl;
        LOG << a*a << nl;
        assert( a*a == 1 + 4 + 9 + 16 + 25 );

        LOG << "Copy the middle slice:" << nl;
        a.setslice( 1, t.getslice(1,3) );
        LOG << a.data_as_text() << nl;
        assert( a[0] == 1 and a[1] == 4 and a[2] == 7 and a[3] == 10 and a[4] == 5 );

        LOG << "Add the middle slice:" << nl;
        a.addslice( 1, t.getslice(1,3), 1000. );
        LOG << a.data_as_text() << nl;
        assert( a[0] == 1 and a[1] == 4004 and a[2] == 7007 and a[3] == 10010 and a[4] == 5 );

        FloatVector e(0);
        LOG << "Should be the zero-dimensional vector:" << nl;
        LOG << e.data_as_text() << nl;
        assert( e.getdimension() == 0 );

        assert( a == a ); assert( not ( a != a ) );
        assert( t != a ); assert( not ( t == a ) );
        assert( e == e ); assert( not ( e != e ) );

    }

    if(true)
    {

        FloatVector a { 1, 3, 0 };
        FloatVector b { -1, -3, 0 };
        FloatVector c { 2, 4, 1 };
        FloatVector d { -5, -4, -3 };
        FloatVector e { 0,0,0 };

        LOG << FloatVector {3} << nl;

        LOG << "positive:     (no ) " << a.is_positive() << nl;
        LOG << "negative:     (no ) " << a.is_negative() << nl;
        LOG << "non-negative: (yes) " << a.is_nonnegative() << nl;
        LOG << "non-positive: (no ) " << a.is_nonpositive() << nl;
        LOG << "zero:         (no ) " << a.is_zero() << nl;

        LOG << "positive:     (no ) " << b.is_positive() << nl;
        LOG << "negative:     (no ) " << b.is_negative() << nl;
        LOG << "non-negative: (no ) " << b.is_nonnegative() << nl;
        LOG << "non-positive: (yes) " << b.is_nonpositive() << nl;
        LOG << "zero:         (no ) " << b.is_zero() << nl;

        LOG << "positive:     (yes) " << c.is_positive() << nl;
        LOG << "negative:     (no ) " << c.is_negative() << nl;
        LOG << "non-positive: (no ) " << c.is_nonpositive() << nl;
        LOG << "non-negative: (yes) " << c.is_nonnegative() << nl;
        LOG << "zero:         (no ) " << c.is_zero() << nl;

        LOG << "positive:     (no ) " << d.is_positive()   << nl;
        LOG << "negative:     (yes) " << d.is_negative() << nl;
        LOG << "non-positive: (yes) " << d.is_nonpositive() << nl;
        LOG << "non-negative: (no ) " << d.is_nonnegative() << nl;
        LOG << "zero:         (no ) " << d.is_zero() << nl;

        LOG << "positive:     (no ) " << e.is_positive()   << nl;
        LOG << "negative:     (no ) " << e.is_negative()   << nl;
        LOG << "non-positive: (yes) " << e.is_nonpositive() << nl;
        LOG << "non-negative: (yes) " << e.is_nonnegative() << nl;
        LOG << "zero:         (yes) " << e.is_zero()       << nl;

        // Assertions for vector a
        assert(a.is_positive() == false);
        assert(a.is_negative() == false);
        assert(a.is_nonnegative() == true);
        assert(a.is_nonpositive() == false);
        assert(a.is_zero() == false);

        // Assertions for vector b
        assert(b.is_positive() == false);
        assert(b.is_negative() == false);
        assert(b.is_nonnegative() == false);
        assert(b.is_nonpositive() == true);
        assert(b.is_zero() == false);

        // Assertions for vector c
        assert(c.is_positive() == true);
        assert(c.is_negative() == false);
        assert(c.is_nonnegative() == true);
        assert(c.is_nonpositive() == false);
        assert(c.is_zero() == false);

        // Assertions for vector d
        assert(d.is_positive() == false);
        assert(d.is_negative() == true);
        assert(d.is_nonnegative() == false);
        assert(d.is_nonpositive() == true);
        assert(d.is_zero() == false);

        // Assertions for vector e
        assert(e.is_positive() == false);
        assert(e.is_negative() == false);
        assert(e.is_nonnegative() == true);
        assert(e.is_nonpositive() == true);
        assert(e.is_zero() == true);

        // LOG << "norm: (3.162)" << a.norm() << space << a.normalize().norm() << nl;
        // LOG << "norm: (3.162)" << b.norm() << space << b.normalize().norm() << nl;
        // LOG << "norm: (4.582)" << c.norm() << space << c.normalize().norm() << nl;
        // LOG << "norm: (7.071)" << d.norm() << space << d.normalize().norm() << nl;
        // LOG << "norm: (0.000)" << e.norm() << nl;

        Float norm_a = a.norm();
        Float norm_b = b.norm();
        Float norm_c = c.norm();
        Float norm_d = d.norm();
        Float norm_e = e.norm();

        // Assertions for norms
        Assert(norm_a == std::sqrt((Float)10), norm_a ); // Approximately std::sqrt(10)
        assert(norm_b == std::sqrt((Float)10) ); // Approximately std::sqrt(10)
        assert(norm_c == std::sqrt((Float)21) ); // Approximately std::sqrt(21)
        assert(norm_d == std::sqrt((Float)50) ); // Approximately std::sqrt(50)
        assert(norm_e == 0. );

        // Normalize vectors
        a.normalize();
        b.normalize();
        c.normalize();
        d.normalize();
        
        // Calculate norms after normalization
        Float norm_a_normalized = a.norm();
        Float norm_b_normalized = b.norm();
        Float norm_c_normalized = c.norm();
        Float norm_d_normalized = d.norm();
        
        // Assertions for norms after normalization
        assert( is_numerically_one( norm_a_normalized ) );
        assert( is_numerically_one( norm_b_normalized ) );
        assert( is_numerically_one( norm_c_normalized ) );
        assert( is_numerically_one( norm_d_normalized ) );
        
    }
    
    {
        // Set entries and check values
        FloatVector w(5,0.);
        LOG << w << nl;
        w.random();
        LOG << w << nl;
    }
        

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
