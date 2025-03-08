
#include <cmath>

#include "../../basic.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test and Benchmark for Factorials and Binomials" << nl;
    
    /* survey factorials */

    LOG << "\nLargest number whose factorial fits into data type (signed)" << nl;
    LOG << "    case signed char       : " << largest_factorial_base<       signed char>() << nl;
    LOG << "    case signed short      : " << largest_factorial_base<      signed short>() << nl;
    LOG << "    case signed int        : " << largest_factorial_base<        signed int>() << nl;
    LOG << "    case signed long       : " << largest_factorial_base<       signed long>() << nl;
    LOG << "    case signed long long  : " << largest_factorial_base<  signed long long>() << nl;

    LOG << "\nLargest number whose factorial fits into data type (unsigned)" << nl;
    LOG << "    case unsigned char     : " << largest_factorial_base<     unsigned char>() << nl;
    LOG << "    case unsigned short    : " << largest_factorial_base<    unsigned short>() << nl;
    LOG << "    case unsigned int      : " << largest_factorial_base<      unsigned int>() << nl;
    LOG << "    case unsigned long     : " << largest_factorial_base<     unsigned long>() << nl;
    LOG << "    case unsigned long long: " << largest_factorial_base<unsigned long long>() << nl;

    /* check factorials */

    LOG << "\nChecking factorial implementations...\n";

    for( int n = 0; n < 50; n++ ) {
        auto x = largest_factorial_base_AUX( n, 2 );
        // printf( "Largest k such that k! <= %d : %lld\n", n, (long long int)x );
        if(  0 <= n && n <=   1 ) Assert( x == 1 );
        if(  2 <= n && n <=   5 ) Assert( x == 2 );
        if(  6 <= n && n <=  23 ) Assert( x == 3 );
        if( 24 <= n && n <= 124 ) Assert( x == 4 );
    }

    LOG << "Integers: ";
    for( int i = 0; i <= largest_factorial_base<int>(); i++ ) {

        int64_t f1 = factorial_integer( i );
        int64_t f2 = factorial_integer_naive( i );
        int64_t f3 = factorial_integer_table( i );
        int64_t f4 = factorial_integer_loop( i );

        LOG << i << space;

        assert( f1 == f2 and f2 == f3 and f3 == f4 );

    }

    LOG << nl;

    LOG << "Float: ";
    for( int i = 0; i < 20; i++ ) {

        Float f1 = factorial_numerical( i );
        Float f2 = factorial_numerical_naive( i );
        Float f3 = factorial_numerical_table( i );
        Float f4 = factorial_numerical_loop( i );

        LOG << i << space;

        assert( is_numerically_close( f1/f2, 1. ) );
        assert( is_numerically_close( f1/f3, 1. ) );
        assert( is_numerically_close( f1/f4, 1. ) );

        assert( is_numerically_close( f1, f2 ) );
        assert( is_numerically_close( f1, f3 ) );
        assert( is_numerically_close( f1, f4 ) );

    }

    LOG << nl;


    LOG << "\nChecking binomial implementations...\n";
    
    for( int n = 0; n <= 7; n++ ) assert( binomial_integer( n, 0 )   == 1 );
    for( int n = 0; n <= 7; n++ ) assert( binomial_integer( n, n )   == 1 );
    for( int n = 1; n <= 7; n++ ) assert( binomial_integer( n, 1 )   == n );
    for( int n = 1; n <= 7; n++ ) assert( binomial_integer( n, n-1 ) == n );

    assert( binomial_integer( 0, 0 ) == 1 );
    
    for( int n = 0; n <= 7; n++ ) 
    for( int k = 1; k <= 6; k++ ) 
    {
        assert( binomial_integer( n, n + k ) == 0 );
    }

    for( int n = 0; n <= 7; n++ ) 
    for( int k = 0; k <= n; k++ ) 
    {
        assert( binomial_integer( n, k ) == factorial_integer(n) / factorial_integer(k) / factorial_integer(n-k) );
    }

    for( int n = 0; n <= 7; n++ ) 
    for( int k = 0; k <= n; k++ ) 
    {
        assert( binomial_integer( n, k ) == binomial_integer_secured( n, k ) );
    }

    assert(  1 == binomial_integer(0, 0) ); 
    assert(  1 == binomial_integer(5, 0) ); 
    assert(  1 == binomial_integer(5, 5) ); 
    assert(  1 == binomial_integer(7, 7) ); 
    assert(  1 == binomial_integer(10, 0)); 
    assert(  6 == binomial_integer(4, 2) ); 
    assert( 10 == binomial_integer(5, 3) ); 
    assert( 10 == binomial_integer(5, 2) ); 
    assert( 15 == binomial_integer(6, 4) ); 
    assert( 20 == binomial_integer(6, 3) ); 
    assert( 70 == binomial_integer(8, 4) ); 
    assert( 84 == binomial_integer(9, 3) ); 
    assert( 35 == binomial_integer(7, 3) ); 
    assert( 28 == binomial_integer(8, 2) ); 
    assert( 45 == binomial_integer(10, 2)); 
    
    assert( 252 == binomial_integer(10, 5) );
    assert( 252 == binomial_integer(10, 5) ); 
    assert( 210 == binomial_integer(10, 4) ); 


    for( int n = 0; n <= 12; n++ )
    for( int k = 0; k <=  n; k++ )
    {
        const auto i = binomial_integer(n,k);
        const auto f = binomial_numerical(n,k);
        Assert( i == f, n, k, i, f );
    }

    for( int n = 1; n <= 12; n++ ) 
    {
        for( int k = 1; k < n; k++ ) 
        {
            // Pascal's identity: C(n, k) = C(n-1, k-1) + C(n-1, k)
            int left  = binomial_integer(n, k);
            int right = binomial_integer(n - 1, k - 1) + binomial_integer(n - 1, k);

            Assert( left == right, left, right );
        }
    }
    
    LOG << "\nAll tests passed successfully!" << nl;
    
    return 0;
}
