#ifndef INCLUDEGUARD_BASIC_HPP
#define INCLUDEGUARD_BASIC_HPP

#include <cassert>     /* assert macro */
#include <cmath>     
#include <cstdint>     
#include <cstdlib>     
#include <ctime>     

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#define unreachable abort 
// __builtin_unreachable

// // #include "assertion.hpp"

#ifndef EXTENDED_PRECISION
typedef double Float;
#else 
typedef long double Float;
#endif

static const Float notanumber = std::numeric_limits<Float>::quiet_NaN();

static const Float machine_epsilon = std::numeric_limits<Float>::epsilon();

static const Float desired_precision = 100 * machine_epsilon;




static const char space = ' ';

static const char* emptystring = "";

static const char nl = '\n';

static const char tab = '\t';






/////////////////////////////////////////////////
//                                             //
//    SIMPLE AUXILIARY ARITHMETICS             //
//                                             //
/////////////////////////////////////////////////



template<typename T>
inline int kronecker( const T& i, const T& j )
{
    if( i == j )
        return 1;
    else
        return 0;
}


template<typename T>
T absolute( const T& n )
{
    assert( n >= 0 or n <= 0 );
    if( n >= 0 )
        return n;
    else
        return -n;
}

template<typename T>
T maximum( const T& a, const T& b )
{
    assert( a >= b or a <= b );
    if( a >= b )
        return a;
    else
        return b;
}

template<typename T>
T minimum( const T& a, const T& b )
{
    assert( a >= b or a <= b );
    if( a <= b )
        return a;
    else
        return b;
}


template<typename T>
T square( const T& x )
{
    return x * x;
}










/////////////////////////////////////////////////
//                                             //
//               POWER FUNCTIONS               //
//                                             //
/////////////////////////////////////////////////


// template<typename T>
// static inline T power( const T& base, const T& exponent )
// {
//     static_assert( not std::is_floating_point<T>::value );
//     assert( base != 0 );
//     assert( exponent >= 0 );
//     if( exponent == 0 ) return 1;
//     return base * power( base, exponent - 1 );
// }

static inline Float power_numerical( Float base, Float exponent )
{
    return std::pow( base, exponent );
}

static inline int power_integer( int base, int exponent )
{
    assert( base != 0 or exponent != 0 );
    assert( exponent >= 0 );
    if( exponent == 0 ) return 1;
    return base * power_integer( base, exponent - 1 );
}

static inline int poweroftwo( int exponent )
{
    return power_integer( 2, exponent );
}

static inline int signpower( int exponent )
{
    return exponent % 2 == 0 ? 1. : -1;
}








/////////////////////////////////////////////////
//                                             //
//        INTEGRAL FACTORIAL, BINOMIALS        //
//               AND AUXILIARIES               //
//                                             //
//   NOTE:                                     //
//   For small inputs, the naive method seems  //
//   to perform best, the table method is only //
//   slightly slower, and the loop method is   //
//   consistently slowest. The differences     //
//   are in the range of 5%, so fairly small   //
//   practically.                              //
//                                             //
/////////////////////////////////////////////////


template<typename T>
inline constexpr uintmax_t largest_factorial_base_AUX( T n, uintmax_t k )
{
    static_assert( std::is_fundamental<T>::value and std::is_integral<T>::value, "T must be a fundamental integral value." );
    
    if( k > n )
        return k-1;
    else
        return largest_factorial_base_AUX<T>( n / T(k), k+1 );
}

template<typename T>
inline constexpr uintmax_t largest_factorial_base()
{
    static_assert( std::is_fundamental<T>::value and std::is_integral<T>::value, "T must be a fundamental integral value." );
    
    const uintmax_t n = std::numeric_limits<T>::max();
    return largest_factorial_base_AUX<T>( n, 2 );
}





static inline uintmax_t factorial_integer_table_old( uintmax_t n )
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

static inline uintmax_t factorial_integer_table( uintmax_t n )
{
    static const uintmax_t facs[21] = {
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

    assert( 0 <= n and n <= 20 );
    return facs[n];
}

static inline uintmax_t factorial_integer_naive( uintmax_t n )
{
    assert( 0 <= n and n <= 20 );
    if( n == 0 ) { 
        return 1;
    } else {
        return n * factorial_integer_naive( n-1 );
    }
}

static inline uintmax_t factorial_integer_loop( uintmax_t n )
{
    assert( 0 <= n and n <= 20 );
    uintmax_t ret = 1;
    while( n > 0 ) ret *= n--;
    return ret;
}






static inline int factorial_integer( int n )
{
    assert( n >= 0 );
    assert( n <= 20 );
    assert( n <= largest_factorial_base<decltype(n)>() );
    
    #ifdef NDEBUG 
    uintmax_t result = factorial_integer_loop( n );
    #else
    uintmax_t result = factorial_integer_table( n );
    #endif
    
    assert( result <= std::numeric_limits<int>::max() );
    return static_cast<int>(result);
}

static inline int binomial_integer( int n, int k )
{
    if( 0 > n ) std::cout << n << std::endl;
    assert( 0 <= n );
    if( k < 0 or n < k )
        return 0;
    uintmax_t result = factorial_integer(n) / ( factorial_integer(k) * factorial_integer(n-k) );
    assert( result <= std::numeric_limits<int>::max() );
    return static_cast<int>(result);
}

















/////////////////////////////////////////////////
//                                             //
//       NUMERICAL FACTORIAL, BINOMIALS        //
//               AND AUXILIARIES               //
//                                             //
//   NOTE:                                     //
//   For small inputs, the loop method seems   //
//   to be fastest, whereas the recursive form //
//   of the faculty performs 10% slower.       //
//                                             //
/////////////////////////////////////////////////

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


static inline Float binomial_numerical( int64_t n, int64_t k )
{
    assert( 0 <= n );
    if( k < 0 or n < k )
        return 0.;
    return factorial_numerical(n) / ( factorial_numerical(k) * factorial_numerical(n-k) );
}














/////////////////////////////////////////////////
//                                             //
//              UNSORTED FUNCTIONS             //
//                                             //
/////////////////////////////////////////////////





template<typename T>
static inline void setmemory( T* pointer, size_t number, const T& value )
{
    assert( pointer != nullptr );
    assert( number >= 0 );
    for( int i = 0; i < number; i++ ) pointer[i] = value;
}



static inline bool issmall( Float value, Float threshold = 100. * std::numeric_limits<Float>::epsilon() )
{
    return absolute(value) < threshold;
}

static inline bool isabout( Float value1, Float value2, Float threshold = 100. * std::numeric_limits<Float>::epsilon() )
{
    return issmall( value1 - value2, threshold );
}






static inline int cast_size_to_int( unsigned long long size )
{
    assert( size < std::numeric_limits<int>::max() );
    return static_cast<int>( size );
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

inline std::vector<int> range( int to )
{
    assert( to >= 0 );
    std::vector<int> ret(to+1);
    for( int i = 0; i <= to; i++ ) ret.at(i) = i;
    assert( ret.size() == to+1 );
    return ret;
}














/////////////////////////////////////////////////
//                                             //
//              STRING OPERATIONS              //
//                                             //
/////////////////////////////////////////////////

inline int count_white_space( const std::string& str ) 
{ 
    int ret = 0;
    
    for( int c = 0; c < str.size(); c++ ) 
        if( isspace( str[c] ) ) 
            ret++; 
    
    return ret;
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









typedef clock_t timestamp;

inline timestamp gettimestamp()
{
    return clock(); 
}

inline std::string timestamp2string( timestamp t )
{
    return std::to_string( static_cast<long double>(t) / CLOCKS_PER_SEC ) + "s";
}









/////////////////////////////////////////////////
//                                             //
//            GAUSSIAN VARIABLES               //
//                                             //
/////////////////////////////////////////////////



// Based on the implementations in the C-FAQ:
// http://c-faq.com/lib/gaussian.html

inline Float gaussrand_1()
{
    const int NSUM = 25;
    
    Float x = 0;
    
    for( int i = 0; i < NSUM; i++) x += rand() / static_cast<Float>(RAND_MAX);
    
    x -= NSUM / 2.0;
    x /= std::sqrt( NSUM / 12.0 );
    
    return x;
}

inline Float gaussrand_2()
{
    static bool phase = false;
    static Float U, V;
    const Float PI = 3.141592654;
    Float Z;

    if( phase ) {
        Z = std::sqrt( -2 * std::log(U) ) * std::cos( 2 * PI * V );
    } else {
        U = ( rand() + 1. ) / ( RAND_MAX + 2. );
        V = rand() / ( RAND_MAX + 1. );
        Z = std::sqrt( -2 * std::log(U) ) * std::sin( 2 * PI * V );
    }
        
    phase = not phase;

    return Z;
}


inline Float gaussrand()
{
    return gaussrand_1();
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

template<typename T>
void mergeelementsinsortedlist
( std::list<T>& L, 
  const std::function<T( const T&, const T& )>& merge,
  const std::function<bool( const T&, const T& )>& compare
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





template <typename T, size_t N>
std::ostream& operator<<( std::ostream& stream, const std::array<T, N>& v)
{
    for( const auto& item : v )
        stream << item << space;
    stream << nl;
    return stream;
}












#endif
