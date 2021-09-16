#ifndef INCLUDEGUARD_BASIC_HPP
#define INCLUDEGUARD_BASIC_HPP


#if __cplusplus < 201703L
#error Compilation of this software requires at least C++14. C++17 is recommended.
#endif


#include <cassert>     /* assert macro */
#include <cmath>     
#include <cstdint>     
#include <cstdio>     
#include <cstdlib>     
#include <ctime>     

#include <algorithm>
#include <array>
#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <list>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>





#include "debug.hpp"



/////////////////////////////////////////////////
//                                             //
//          FLOATING POINT DEFINITIONS         //
//                                             //
/////////////////////////////////////////////////

#ifndef EXTENDED_PRECISION
typedef double Float;
#else 
typedef long double Float;
#endif

inline const constexpr Float notanumber = std::numeric_limits<Float>::quiet_NaN();

inline const constexpr Float machine_epsilon = std::numeric_limits<Float>::epsilon();

inline const constexpr Float desired_precision = 100. * machine_epsilon;








/////////////////////////////////////////////////
//                                             //
//        CHAR AND STRING CONSTANTS            //
//                                             //
/////////////////////////////////////////////////

inline const constexpr char space = ' ';

inline const constexpr char* emptystring = "";

inline const constexpr char nl = '\n';

inline const constexpr char tab = '\t';








/////////////////////////////////////////////////
//                                             //
//        SIMPLE AUXILIARY ARITHMETICS         //
//                                             //
/////////////////////////////////////////////////



template<typename T>
inline constexpr int kronecker( const T& i, const T& j )
{
    if( i == j )
        return 1;
    else
        return 0;
}


template<typename T>
inline constexpr T absolute( const T& n )
{
    if( n >= 0 ) return  n;
    if( n <= 0 ) return -n;
    assert( not std::isfinite(n) );
    return n;
}

template<typename T>
inline constexpr T maximum( const T& a, const T& b )
{
    if( a >= b ) return a;
    if( a <= b ) return b;
    assert( ( not std::isfinite(a) ) or ( not std::isfinite(b) ) );
    return a;
//     assert( a >= b or a <= b );
//     if( a >= b )
//         return a;
//     else
//         return b;
}

template<typename T>
inline constexpr T minimum( const T& a, const T& b )
{
    if( a <= b ) return a;
    if( a >= b ) return b;
    assert( ( not std::isfinite(a) ) or ( not std::isfinite(b) ) );
    return a;
//     assert( a >= b or a <= b );
//     if( a >= b )
//         return a;
//     else
//         return b;
}


template<typename T>
inline constexpr T square( const T& x )
{
    return x * x;
}



inline constexpr bool issmall( Float value, Float threshold = 100. * machine_epsilon )
{
    return absolute(value) < threshold;
}

inline constexpr bool isaboutequal( Float value1, Float value2, Float threshold = 100. * machine_epsilon )
{
    return issmall( value1 - value2, threshold );
}








/////////////////////////////////////////////////
//                                             //
//               POWER FUNCTIONS               //
//                                             //
/////////////////////////////////////////////////


// template<typename T>
// inline constexpr T power( const T& base, const T& exponent )
// {
//     static_assert( not std::is_floating_point<T>::value );
//     assert( base != 0 );
//     assert( exponent >= 0 );
//     if( exponent == 0 ) return 1;
//     return base * power( base, exponent - 1 );
// }

inline /*constexpr*/ Float power_numerical( Float base, Float exponent )
{
    return std::pow( base, exponent );
}

inline constexpr int power_integer( int base, int exponent )
{
    assert( base != 0 or exponent != 0 );
    assert( exponent >= 0 );
    if( exponent == 0 ) return 1;
    return base * power_integer( base, exponent - 1 );
}

inline constexpr int poweroftwo( int exponent )
{
    return power_integer( 2, exponent );
}

inline constexpr int signpower( int exponent )
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

/*
 * Recursively divide the integer n by larger and larger numbers 1, 2, 3, ... without remainder
 * until the divisor is larger than n. That divisor is the largest numbers 
 * whose factorial is at most n.
 */

template<typename T>
inline constexpr uintmax_t largest_factorial_base_AUX( T n, uintmax_t k )
{
    static_assert( std::is_fundamental<T>::value and std::is_integral<T>::value, "T must be a fundamental integral value." );

#if __cplusplus >= 201402L
    if( k > n )
        return k-1;
    else
        return largest_factorial_base_AUX<T>( n / T(k), k+1 );
#else
    return (k > n) ? (k-1) : (largest_factorial_base_AUX<T>( n / T(k), k+1 ));
#endif
}

template<typename T>
inline constexpr uintmax_t largest_factorial_base()
{
    static_assert( std::is_fundamental<T>::value and std::is_integral<T>::value, "T must be a fundamental integral value." );
    
#if __cplusplus >= 201402L
    const uintmax_t n = std::numeric_limits<T>::max();
    return largest_factorial_base_AUX<T>( n, 2 );
#else
    return largest_factorial_base_AUX<T>( std::numeric_limits<T>::max(), 2 );
#endif
}





inline constexpr uintmax_t factorial_integer_table_old( uintmax_t n )
{
    switch(n){
        case 0: 
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

inline constexpr uintmax_t factorial_integer_table( uintmax_t n )
{
    constexpr uintmax_t facs[21] = {
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

inline constexpr uintmax_t factorial_integer_naive( uintmax_t n )
{
    assert( 0 <= n and n <= 20 );
    if( n == 0 ) { 
        return 1;
    } else {
        return n * factorial_integer_naive( n-1 );
    }
}

inline constexpr uintmax_t factorial_integer_loop( uintmax_t n )
{
    assert( 0 <= n and n <= 20 );
    uintmax_t ret = 1;
    while( n > 0 ) ret *= n--;
    return ret;
}






inline constexpr int factorial_integer( int n )
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

inline constexpr int binomial_integer( int n, int k )
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

inline constexpr Float factorial_numerical_naive( int64_t n )
{
    assert( 0 <= n );
    if( n == 0 ) { 
        return 1.;
    } else {
        return static_cast<Float>(n) * factorial_numerical_naive( n-1 );
    }
}

inline constexpr Float factorial_numerical_loop( int64_t n )
{
    assert( 0 <= n );
    Float ret = 1.;
    while( n > 0 ) ret *= static_cast<Float>(n--);
    return ret;
}

inline constexpr Float factorial_numerical_table( int64_t n )
{
    constexpr Float facs[21] = {
        1.,
        1.,
        2.,
        6.,
        24.,
        120.,
        720.,
        5040.,
        40320.,
        362880.,
        3628800,
        39916800.,
        479001600.,
        6227020800.,
        87178291200.,
        1307674368000.,
        20922789888000.,
        355687428096000.,
        6402373705728000.,
        121645100408832000.,
        2432902008176640000.,
    };

    assert( 0 <= n and n <= 20 );
    return facs[n];
}

inline constexpr Float factorial_numerical( int64_t n )
{
    return factorial_numerical_table( n );
}


inline constexpr Float binomial_numerical( int64_t n, int64_t k )
{
    assert( 0 <= n );
    if( k < 0 or n < k )
        return 0.;
    return factorial_numerical(n) / ( factorial_numerical(k) * factorial_numerical(n-k) );
}










/////////////////////////////////////////////////
//                                             //
//            GAUSSIAN VARIABLES               //
//                                             //
/////////////////////////////////////////////////




inline void seed_random_integer()
{
    srand(0);
}

inline int random_integer()
{
    int ret = rand();
    assert( 0 <= ret and ret <= RAND_MAX );
    return ret;
}

inline Float random_uniform()
{
    Float ret = rand() / static_cast<Float>( RAND_MAX );
    assert( 0. <= ret and ret <= 1. );
    return ret;
}


// Based on the implementations in the C-FAQ:
// http://c-faq.com/lib/gaussian.html

inline Float gaussrand_1()
{
    const int NSUM = 25;
    
    Float x = 0;
    
    for( int i = 0; i < NSUM; i++) 
        x += rand() / static_cast<Float>(RAND_MAX);
    
    x -= NSUM / 2.0;
    x /= std::sqrt( NSUM / 12.0 );
    
    return x;
}

inline Float gaussrand_2()
{
    static bool phase = false;
    static Float U, V;
    const Float PI = 3.14159265358979323846;
    Float Z;

    if( phase ) {
        Z = std::sqrt( -2. * std::log(U) ) * std::cos( 2. * PI * V );
    } else {
        U = ( rand() + 1. ) / ( RAND_MAX + 2. );
        V = rand() / ( RAND_MAX + 1. );
        Z = std::sqrt( -2. * std::log(U) ) * std::sin( 2. * PI * V );
    }
        
    phase = not phase;

    return Z;
}


// http://c-faq.com/lib/gaussrand.luben.html
inline Float gaussrand_3( Float mean = 0., Float std_dev = 1. )
{
    assert( std_dev > machine_epsilon );

    Float x = rand() / (RAND_MAX + 1.0);   /* 0.0 <= y < 1.0 */
    
    unsigned low = (x < 0.5) ? 0 : 1;
    
    Float y = std::abs(x - 1.0);                        /* 0.0 < y <= 1.0 */
    Float z = std_dev * std::sqrt( -2.0 * std::log(y) );

    return low ? (mean + z) : (mean - z);
}


inline Float gaussrand()
{
    return gaussrand_1();
}





/*
inline void random_unit_vector( Float* values, const int N )
{
    assert( N >= 0 );
    assert( values != nullptr );
    
    const Float PI = 3.14159265358979323846;
    
    Float norm_sq = 0;
    
    for( int k = 0; k < N/2; k++ ) {
        
        int k1 = 2*k;
        int k2 = k1+1;
        
        Float U = ( rand() + 1. ) / ( RAND_MAX + 2. );
        Float V = rand() / ( RAND_MAX + 1. );
        
        Float radius_sq = -2. * std::log( U );
        Float radius = std::sqrt( radius_sq );
        
        values[k1] = radius * std::sin( 2. * PI * V );
        values[k2] = radius * std::cos( 2. * PI * V );
        
        norm_sq += radius_sq;
    }
    
    if( N % 2 == 1 ) {
        values[N-1] = gaussrand();
        norm_sq += values[N-1] * values[N-1];
    }
    
    
}
*/















/////////////////////////////////////////////////
//                                             //
//              TIME UTILITIES                 //
//                                             //
/////////////////////////////////////////////////


static_assert( std::is_integral< decltype( std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count() ) >::value , "Time measurement must be integral" );


typedef uintmax_t timestamp;

inline timestamp gettimestamp()
{
    
    static timestamp start_timestamp = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    timestamp now = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    assert( now >= start_timestamp );
    
    return now - start_timestamp;
}



inline std::string timestamp2measurement( const timestamp& t )
{
    return std::to_string( static_cast<long double>(t) ) + "ms";
}

// inline std::string measurementnow( const timestamp& t ) // TODO Remove this line 
inline std::string measurementnow()
{
    return timestamp2measurement( gettimestamp() );
}



inline std::string timestamp2digitalcode( const timestamp& t )
{
    std::ostringstream ss;
//     ss.reserve(12);
//     ss << std::setw(14) << t;
    ss << std::hex << std::setfill('_') << std::setw(8) << t;
    return ss.str();
}

inline std::string digitalcodenow()
{
    return timestamp2digitalcode( gettimestamp() );
}



inline std::string protocolprefixnow()
{
    // static const std::string foo = std::string("\e[36m[");
    // static const std::string bar = std::string("]\e[39m\t");
    static const std::string foo = std::string("[");
    static const std::string bar = std::string("]\t");
    return foo + digitalcodenow() + bar;
}

















/////////////////////////////////////////////////
//                                             //
//       SUM INTEGERS PRODUCED BY LAMBDA       //
//                                             //
/////////////////////////////////////////////////

inline int sum_int( int from, int to, const std::function<int(int)>& calc )
{
    if( from > to )
        return 0;
    int ret = 0;
    for( int i = from; i <= to; i++ )
        ret += calc( i );
    return ret;
}

inline int sum_int( int to, const std::function<int(int)>& calc )
{
    return sum_int( 0, to, calc );
}

























/////////////////////////////////////////////////
//                                             //
//              BUMP FUNCTIONS                 //
//                                             //
/////////////////////////////////////////////////






inline Float bumpfunction( Float x )
{
    Float delta = x*x - 1.;

    if( absolute(x) < 0.99999999 ) {

        return std::exp( 1. / delta );

    } else {

        return 0;
        
    }
}

inline Float bumpfunction_dev( Float x )
{
    
    Float delta = x*x - 1.;
    Float delta_sq = delta*delta;

    if( absolute(x) < 0.99999999 ) {
        
        return -2. * x * std::exp( 1. / delta ) / delta_sq;
        
    } else {
        
        return 0;
        
    }
}

inline Float bumpfunction_devdev( Float x )
{
    
    
//     Float t1 = std::exp( -1. / ( 1. - x*x ) );
//     Float t2 = std::exp( 1 - x*x );
// 
//     if( x*x < 1 )
//         return
//             t1 * ( 4*x*x*pow(t2,-4.) - 2*pow(t2,-2.) - 8*x*x*pow(t2,-3.) );
//     else
//         return
//             0.;
                            

    
    
    Float delta = x*x - 1.;
    
    Float delta_sq = delta    * delta;
    Float delta_p4 = delta_sq * delta_sq;

    if( absolute(x) < 0.99999999 ) {
        
        return std::exp( 1. / delta ) * (  6.*x*x*x*x - 2. ) / delta_p4;
        
    } else {
        
        return 0.;
        
    }
}





















/////////////////////////////////////////////////
//                                             //
//       CARTESIAN AND POLAR COORDINATES       //
//                                             //
/////////////////////////////////////////////////

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














































/////////////////////////////////////////////////
//                                             //
//          MISC LIBRARY HACKS                 //
//                                             //
/////////////////////////////////////////////////



/******************************************************/
/*    use this to safely cast size_types to C++ int   */
/******************************************************/

inline int SIZECAST( std::uintmax_t size )
{
    assert( size < std::numeric_limits<int>::max() );
    return static_cast<int>( size );
}


/******************************************************/
/*      insert tabs before each line       */
/******************************************************/

inline std::string tab_each_line( std::string str ) 
{ 
    str.insert( 0, 1, '\t' );
    for( int c = str.size(); c > 0; c-- ) {
        if( str[c-1] == '\n' )
            str.insert(c, 1, '\t');
    }
    return str;
} 


/******************************************************/
/*      count the white space within STL string       */
/******************************************************/

inline int count_white_space( const std::string& str ) 
{ 
    int ret = 0;
    
    for( int c = 0; c < str.size(); c++ ) 
        if( isspace( str[c] ) ) 
            ret++; 
    
    return ret;
} 


/******************************************************/
/*        write zero-based range into vector          */
/******************************************************/

inline std::vector<int> range( int to )
{
    assert( to >= 0 );
    std::vector<int> ret(to+1);
    for( int i = 0; i <= to; i++ ) ret.at(i) = i;
    assert( ret.size() == to+1 );
    return ret;
}


/******************************************************/
/*   remove duplicates from random access container   */
/******************************************************/

template< typename T >
inline void sort_and_remove_duplicates( T& t )
{
    std::sort( t.begin(), t.end() );
    auto last = std::unique( t.begin(), t.end() );
    t.erase( last, t.end() );
}


/******************************************************/
/*      find index of element with STL vector         */
/******************************************************/

template<typename T>
int find_index( const std::vector<T>& vec, const T& t )
{
   const auto it = std::find( vec.begin(), vec.end(), t );
   assert( it != vec.end() );
   const auto ret = std::distance( vec.begin(), it );
   assert( ret >= 0 );
   assert( ret < vec.size() );
   return SIZECAST( ret );
}


/******************************************************/
/*            merge two sorted STL lists              */
/******************************************************/

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




/***********************************************/
/*   GENERIC STREAM TEMPLATE FOR std::array    */ 
/***********************************************/

template <typename T, size_t N>
std::ostream& operator<<( std::ostream& stream, const std::array<T, N>& v)
{
    for( const auto& item : v )
        stream << item << space;
    stream << nl;
    return stream;
}




/***********************************************/
/*         make_unique HACK                    */ 
/***********************************************/

#if __cplusplus < 201402L

/****
 * 
 * A very imperfect solution for make_unique in C++11
 * We enter undefined behavior territory here
 * 
 ****/
#warning \
This code extends the std namespace so that `make_unique` \
is available throughout the code. This was triggered by a C++ version \
below C++14. While this may be a practical workaround, it is officially \
considered undefined behavior in the C++ standard. Please try to compile \
with C++14 or higher.

#include <memory>

namespace std
{
template <typename T, typename ...Args> 
std::unique_ptr<T> make_unique(Args && ...args)
{
  return std::unique_ptr<T>( new T(std::forward<Args>(args)...) );
}
}

#endif









// template<typename T>
// inline void setmemory( T* pointer, size_t number, const T& value )
// {
//     assert( pointer != nullptr );
//     assert( number >= 0 );
//     for( int i = 0; i < number; i++ ) pointer[i] = value;
// }



// inline void sort_integers( int* start, int length )
// {
//     assert( start != nullptr && length >= 0 );
//     for( int i = 1; i < length; i++ )
//     for( int j = 1; j < length; j++ )
//         if( start[j-1] > start[j] ) 
//             std::swap( start[j-1], start[j] );
// }









#endif
