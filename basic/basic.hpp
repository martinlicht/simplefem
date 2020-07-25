#ifndef INCLUDEGUARD_BASIC_HPP
#define INCLUDEGUARD_BASIC_HPP

#include <cstdint>     
#include <cmath>     
#include <ctime>     
#include <cstdlib>     
#include <cassert>     /* assert macro */
#include <list>
#include <iterator>
#include <functional>
#include <iostream>

#define unreachable abort 
// __builtin_unreachable

// // #include "assertion.hpp"

typedef long double Float;

static const Float notanumber = std::numeric_limits<Float>::quiet_NaN();




static const char space = ' ';

static const char* emptystring = "";

static const char nl = '\n';

static const char tab = '\t';





template<typename T>
inline int kronecker( T i, T j )
{
    if( i == j )
        return 1;
    else
        return 0;
}


template<typename T>
T absolute( const T& n )
{
    if( n >= 0 )
        return n;
    else
        return -n;
}

template<typename T>
T maximum( const T& a, const T& b )
{
    if( a >= b )
        return a;
    else
        return b;
}

template<typename T>
T minimum( const T& a, const T& b )
{
    if( a <= b )
        return a;
    else
        return b;
}






static inline int64_t factorial_integer_table_old( int64_t n )
{
    switch(n){
        case 0: return 1;
        case 1: return 1;
        case 2: return 2;
        case 3: return 6;
        case 4: return 24;
        case 5: return 120;
        case 6: return 720ll;
        case 7: return 5040ll;
        case 8: return 40320ll;
        case 9: return 362880ll;
        case 10: return 3628800ll;
        case 11: return 39916800ll;
        case 12: return 479001600ll;
        case 13: return 6227020800ll;
        case 14: return 87178291200ll;
        case 15: return 1307674368000ll;
        case 16: return 20922789888000ll;
        case 17: return 355687428096000ll;
        case 18: return 6402373705728000ll;
        case 19: return 121645100408832000ll;
        case 20: return 2432902008176640000ll;
        // case 21: return 51090942171709440000ll;
        // case 22: return 1124000727777607680000ll;
        // case 23: return 25852016738884976640000ll;
        // case 24: return 620448401733239439360000ll;
        // case 25: return 15511210043330985984000000ll;
        default: unreachable();
    }
}

static inline int64_t factorial_integer_table( int64_t n )
{
    static const int64_t facs[21] = {
        1,
        1,
        2,
        6,
        24,
        120,
        720ll,
        5040ll,
        40320ll,
        362880ll,
        3628800ll,
        39916800ll,
        479001600ll,
        6227020800ll,
        87178291200ll,
        1307674368000ll,
        20922789888000ll,
        355687428096000ll,
        6402373705728000ll,
        121645100408832000ll,
        2432902008176640000ll,
    };

    assert( 0 <= n and n <= 21 );
    return facs[n];
}

static inline int64_t factorial_integer_naive( int64_t n )
{
    assert( 0 <= n and n <= 21 );
    if( n == 0 ) { 
        return 1;
    } else {
        return n * factorial_integer_naive( n-1 );
    }
}

static inline int64_t factorial_integer_loop( int64_t n )
{
    assert( 0 <= n and n <= 21 );
    int64_t ret = 1;
    while( n > 0 ) ret *= n--;
    return ret;
}

static inline int64_t factorial_integer( int64_t n )
{
    #ifdef NDEBUG 
    return factorial_integer_loop( n );
    #else
    return factorial_integer_table( n );
    #endif
}







static inline Float factorial_numerical_naive( int64_t n )
{
    assert( 0 <= n );
    if( n == 0 ) { 
        return 1.;
    } else {
        return n * factorial_numerical_naive( n-1 );
    }
}

static inline Float factorial_numerical_loop( int64_t n )
{
    assert( 0 <= n );
    Float ret = 1.;
    while( n > 0 ) ret *= n--;
    return ret;
}

static inline Float factorial_numerical( int64_t n )
{
    return factorial_numerical_loop( n );
}







static inline int64_t binomial_integer( int64_t n, int64_t k )
{
    assert( 0 <= n );
//     assert( 0 <= k && k <= n );
    if( k < 0 or n < k )
        return 0;
    return factorial_integer(n) / ( factorial_integer(k) * factorial_integer(n-k) );
}

static inline Float binomial_numerical( int64_t n, int64_t k )
{
    assert( 0 <= n );
//     assert( 0 <= k && k <= n );
    if( k < 0 or n < k )
        return 0.;
    return factorial_numerical(n) / ( factorial_numerical(k) * factorial_numerical(n-k) );
}



// todo: deprecate these two templated functions

// template<typename T>
// T factorial( const T& n )
// {
//     if( n == 0 )
//         return 1;
//     else if( n < 0 )
//         { unreachable(); }
//     else
//         return n * factorial<T>(n-1);
// }

// template<typename T>
// T binomial( const T& n, const T& k )
// {
//     assert( 0 <= n );
//     assert( 0 <= k && k <= n );
//     return factorial(n) / ( factorial(k) * factorial(n-k) );
// }


static inline Float power( Float base, Float exponent )
{
    return std::pow( base, exponent );
}

template<typename T>
static inline T power( T base, T exponent )
{
    assert( base != 0 );
    assert( exponent >= 0 );
    if( exponent == 0 ) return 1;
    return base * power( base, exponent - 1 );
}



static inline int integerpower( int base, int exponent )
{
    assert( base != 0 );
    assert( exponent >= 0 );
    if( exponent == 0 ) return 1;
    return base * integerpower( base, exponent - 1 );
}

static inline int poweroftwo( int exponent )
{
    return integerpower( 2, exponent );
}

static inline int signpower( int exponent )
{
    return exponent % 2 == 0 ? 1. : -1;
}





static inline int getbit( unsigned int value, unsigned int bitnumber )
{
    return ( value >> bitnumber ) % 2;
}






static inline bool issmall( Float value, Float threshold = 0.00001 )
{
    return absolute(value) < threshold;
}

static inline bool isabout( Float value1, Float value2, Float threshold = 0.00001 )
{
    return issmall( value1 - value2, threshold );
}






static inline int sum_int( int from, int to, const std::function<int(int)>& calc )
{
    if( from > to )
        return 0;
    int ret = 0;
    for( int i = from; i <= to; i++ )
        ret += calc( i );
    return ret;
}

static inline int sum_int( int to, const std::function<int(int)>& calc )
{
    return sum_int( 0, to, calc );
}







typedef clock_t timestamp;

inline timestamp gettimestamp()
{
    return CLOCKS_PER_SEC * static_cast<double>(clock()); 
}


template<typename T>
void mergeelementsinsortedlist
( std::list<T>& L, 
  std::function<T( const T&, const T& )> merge,
  std::function<bool( const T&, const T& )> compare
) {
    typename std::list<T>::iterator it = L.begin();
    while( it != L.end() ){

        typename std::list<T>::iterator now = it; 
        typename std::list<T>::iterator next = ++it;

        if( next == L.end() ) return;

        if( compare( *it, *next ) )
        {
            *now = merge( *now, *next );
            L.erase( next );
            it = now;
        } 

    }
}



#include <vector>
#include <algorithm>

template<typename T>
int find_index( const std::vector<T>& vec, const T& t )
{
   const auto& it = std::find( vec.begin(), vec.end(), t );
   assert( it != vec.end() );
   int ret = std::distance( vec.begin(), it );
   assert( ret >= 0 );
   assert( ret < vec.size() );
   return ret;
}


#include <ostream>
#include <array>


template <typename T, size_t N>
std::ostream& operator<<( std::ostream& stream, const std::array<T, N>& v)
{
    for( const auto& item : v )
        stream << item << space;
    stream << nl;
    return stream;
}




inline Float gaussrand()
{
    const int NSUM = 25;
    
    Float x = 0;
    
    for( int i = 0; i < NSUM; i++) x += rand() / (Float)RAND_MAX;
    
    x -= NSUM / 2.0;
    x /= sqrt( NSUM / 12.0 );
    
    return x;
}



inline void cartesian_to_polar_coordinates2D( const Float& x, const Float& y, Float& radius, Float& angle )
{
    radius = std::sqrt( x*x + y*y );
    angle  = std::atan2( x, y );
}

inline void polar_to_cartesian_coordinates2D( const Float& radius, const Float& angle, Float& x, Float& y )
{
    x = radius * std::cos( angle );
    y = radius * std::sin( angle );
}






inline void sort_integers( int* start, int length )
{
    assert( start != nullptr && length >= 0 );
    for( int i = 1; i < length; i++ )
    for( int j = 1; j < length; j++ )
        if( start[j-1] > start[j] ) 
            std::swap( start[j-1], start[j] );
}



template< typename T >
inline void sort_and_unique( T& t )
{
    std::sort( t.begin(), t.end() );
    auto last = std::unique( t.begin(), t.end() );
    t.erase( last, t.end() );
}






inline std::vector<int> range( int to )
{
    assert( to >= 0 );
    std::vector<int> ret(to+1);
    for( int i = 0; i <= to; i++ ) ret.at(i) = i;
    assert( ret.size() == to+1 );
    return ret;
}







#endif
