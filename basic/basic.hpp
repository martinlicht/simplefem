#ifndef INCLUDEGUARD_BASIC_HPP
#define INCLUDEGUARD_BASIC_HPP


#if __cplusplus < 201402L
#error Compilation of this software requires at least C++14.
#endif


#include <cmath>
// #include <cstdint>
// #include <cstdio>
// #include <cstdlib>
// #include <ctime>

// #include <algorithm>
#include <array>
// #include <chrono>
// #include <functional>
// #include <ostream>
// #include <iterator>
#include <limits>
// #include <list>
#include <string>
#include <type_traits>
// #include <utility>
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


// inline constexpr Float HeronStep( Float x, Float a )
// {
//     return 0.5 * ( x + a / x );
// }

// inline constexpr Float HeronAlgorithm( Float x, Float a, unsigned int i )
// {
//     return i==0 ? x : HeronStep( HeronAlgorithm(x,a,i-1), a );
// }

// template<unsigned int i = 10>
// static Float Sqrt( Float a )
// {
//     return HeronAlgorithm( x, a, i );
// }


static const constexpr Float notanumber = std::numeric_limits<Float>::quiet_NaN();

static const constexpr Float machine_epsilon = std::numeric_limits<Float>::epsilon();

static const /*constexpr*/ Float desired_precision = 
                                    sizeof(Float) == sizeof(float) ? 10e-6 : std::sqrt( machine_epsilon );


// TODO: Put these somewhere where it makes sense 
extern template class std::vector<int>;
extern template class std::vector<Float>;




/////////////////////////////////////////////////
//                                             //
//        CHAR AND STRING CONSTANTS            //
//                                             //
/////////////////////////////////////////////////

static const constexpr char space = ' ';

static const constexpr char* emptystring = "";

static const constexpr char nl = '\n';

static const constexpr char tab = '\t';








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
inline constexpr T absolute( const T& x )
{
    if( x >= 0 ) return  x;
    if( x <= 0 ) return -x;
    Assert( not std::isfinite(x) );
    return x;
}

template<typename T>
inline constexpr T sign( const T& x )
{
    Assert( std::isfinite(x) );
    if( x > 0 ) return  1;
    if( x < 0 ) return -1;
    else        return  0;
}

template<typename T>
inline constexpr T maximum( const T& a, const T& b )
{
    if( a >= b ) return a;
    if( a <= b ) return b;
    Assert( ( not std::isfinite(a) ) or ( not std::isfinite(b) ) );
    return a;
//     Assert( a >= b or a <= b ); if( a >= b ) return a; else return b;
}

template<typename T>
inline constexpr T maxabs( const T& a, const T& b )
{
    return maximum( absolute(a), absolute(b) );
}

template<typename T>
inline constexpr T minimum( const T& a, const T& b )
{
    if( a <= b ) return a;
    if( a >= b ) return b;
    Assert( ( not std::isfinite(a) ) or ( not std::isfinite(b) ) );
    return a;
//     Assert( a >= b or a <= b ); if( a >= b ) return b; else return a;
}


template<typename T, typename... Args >
inline constexpr T maximum( T t, Args... args )
{
    return maximum( t, maximum( args... ) );
}

template<typename T, typename... Args >
inline constexpr T maxabs( T t, Args... args )
{
    return maxabs( t, maxabs( args... ) );
}

template<typename T, typename... Args >
inline constexpr T minimum( T t, Args... args )
{
    return minimum( t, minimum( args... ) );
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

// https://codingnest.com/the-little-things-comparing-floating-point-numbers/






/////////////////////////////////////////////////
//                                             //
//               POWER FUNCTIONS               //
//                                             //
/////////////////////////////////////////////////


inline /*constexpr*/ Float power_numerical( Float base, Float exponent )
{
    return std::pow( base, exponent );
}

inline constexpr int power_integer( int base, int exponent )
{
    Assert( base != 0 or exponent != 0 );
    Assert( exponent >= 0 );
    if( exponent == 0 ) return 1;
    if( base == 0 ) return 0;
    int rec = power_integer( base, exponent - 1 );
    int ret = base * rec;
    Assert( ret / base == rec );
    return ret;
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

    Assert( 0 <= n and n <= 20 );
    return facs[n];
}

inline constexpr uintmax_t factorial_integer_naive( uintmax_t n )
{
    Assert( 0 <= n and n <= 20 );
    if( n == 0 ) { 
        return 1;
    } else {
        return n * factorial_integer_naive( n-1 );
    }
}

inline constexpr uintmax_t factorial_integer_loop( uintmax_t n )
{
    Assert( 0 <= n and n <= 20 );
    uintmax_t ret = 1;
    while( n > 0 ) ret *= n--;
    return ret;
}






inline constexpr int factorial_integer( int n )
{
    Assert( n >= 0 );
    Assert( n <= 20 );
    Assert( n <= largest_factorial_base<decltype(n)>() );
    
    #ifdef NDEBUG 
    uintmax_t result = factorial_integer_table( n );
    #else
    uintmax_t result = factorial_integer_table( n );
    #endif
    
    Assert( result <= std::numeric_limits<int>::max() );
    return static_cast<int>(result);
}

inline constexpr int binomial_integer( int n, int k )
{
    Assert( 0 <= n, "Negative n for integer binomial: ", n ); 
    Assert( 0 <= n );
    if( k < 0 or n < k )
        return 0;
    uintmax_t result = factorial_integer(n) / ( factorial_integer(k) * factorial_integer(n-k) );
    Assert( result <= std::numeric_limits<int>::max() );
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
    Assert( 0 <= n );
    if( n == 0 ) { 
        return 1.;
    } else {
        return static_cast<Float>(n) * factorial_numerical_naive( n-1 );
    }
}

inline constexpr Float factorial_numerical_loop( int64_t n )
{
    Assert( 0 <= n );
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

    Assert( 0 <= n and n <= 20 );
    return facs[n];
}

inline constexpr Float factorial_numerical( int64_t n )
{
    return factorial_numerical_naive( n );
}


inline constexpr Float binomial_numerical( int64_t n, int64_t k )
{
    Assert( 0 <= n );
    if( k < 0 or n < k )
        return 0.;
    return factorial_numerical(n) / ( factorial_numerical(k) * factorial_numerical(n-k) );
}



























/////////////////////////////////////////////////
//                                             //
//              TIME UTILITIES                 //
//                                             //
/////////////////////////////////////////////////



typedef uintmax_t timestamp;

timestamp gettimestamp();

// TODO: move to utility 

std::string timestamp2measurement( const timestamp& t );

// inline std::string measurementnow( const timestamp& t ) // TODO Remove this line 
std::string measurementnow();

std::string timestamp2digitalcode( const timestamp& t );

std::string digitalcodenow();

std::string protocolprefixnow();









/////////////////////////////////////////////////
//                                             //
//              BUMP FUNCTIONS                 //
//                                             //
/////////////////////////////////////////////////

Float bumpfunction( Float x );

Float bumpfunction_dev( Float x );

Float bumpfunction_devdev( Float x );

Float bumpfunction_devdevdev( Float x );

Float blob( Float x );

Float blob_dev( Float x );

Float blob_devdev( Float x );

Float blob_devdevdev( Float x );

Float sinpy( Float x );

Float cospy( Float x );






/////////////////////////////////////////////////
//                                             //
//       CARTESIAN AND POLAR COORDINATES       //
//                                             //
/////////////////////////////////////////////////

void cartesian_to_polar_coordinates2D( const Float& x, const Float& y, Float& radius, Float& angle );

void polar_to_cartesian_coordinates2D( const Float& radius, const Float& angle, Float& x, Float& y );














/////////////////////////////////////////////////
//                                             //
//            STRING UTILITIES                 //
//                                             //
/////////////////////////////////////////////////

/******************************************************/
/*      insert tabs before each line                  */
/*      pad each line                                 */
/******************************************************/

inline std::string tab_each_line( std::string str ) // TODO: Move to utilities 
{ 
    str.insert( 0, 1, '\t' );
    for( int c = str.size(); c > 0; c-- ) {
        if( str[c-1] == '\n' )
            str.insert(c, 1, '\t');
    }
    return str;
} 

// inline std::string pad_each_line( std::string str, std::string pad )
// {
//     int c = 0;
//     for( int i = 0; i < str.length(); i++ ) if( str[i] == '\n' ) c++;
//     std::string ret = pad;
//     ret.reserve( str.length() + c * pad.length() );
//     for( int i = 0; i < str.length(); i++ )
//     {
//         ret += str[i];
//         if( str[i] == '\n' ) ret += pad.length();
    
//     return ret;
// }


/******************************************************/
/*      count the white space within STL string       */
/******************************************************/

inline int count_white_space( const std::string& str ) // TODO: Move to utilities 
{ 
    int ret = 0;  
    for( int c = 0; c < str.size(); c++ ) if( isspace( str[c] ) ) ret++; 
    return ret;
} 











/////////////////////////////////////////////////
//                                             //
//            GENERIC STREAMING                //
//                                             //
/////////////////////////////////////////////////


/***********************************************/
/*   GENERIC STREAM TEMPLATE FOR std::array    */ 
/***********************************************/

// TODO: Move into separate include file 
template <typename S, typename T, size_t N>
inline S& operator<<( S& stream, const std::array<T, N>& v)
{
    for( const auto& item : v )
        stream << item << space;
    stream << nl;
    return stream;
}

/******************************************************/
/*       printf into C++ strings and streams          */
/******************************************************/

// template< typename... Params >
// inline std::string printf_into_string( const char* formatstring, Params... args )
// {
//     std::size_t length = std::snprintf(nullptr, 0, formatstring, args... ) + 1;
//     char* c_str = new char[length];
//     std::snprintf( c_str, length, formatstring, args... );
//     std::string ret( c_str );
//     delete[] c_str;
//     return ret;
// }

std::string printf_into_string( const char* formatstring, ... ) 
__attribute__ (( format (printf,1,2) ));



// template< typename L, typename... Params >
// inline void printf_into_stream( L& stream, const char* formatstring, Params... args )
// {
//     stream << printf_into_string( formatstring, args... );
// }
#if __cplusplus < 202002L
#define printf_into_stream( stream, ... ) { stream << printf_into_string( __VA_ARGS__ ); }
#else
#define printf_into_stream( stream, formatstring, ... ) { stream << printf_into_string( formatstring __VA_OPT__(,) __VA_ARGS__ ); }
#endif



















































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
    Assert( size < std::numeric_limits<int>::max() );
    return static_cast<int>( size );
}





// template<typename T>
// inline void setmemory( T* pointer, size_t number, const T& value )
// {
//     Assert( pointer != nullptr );
//     Assert( number >= 0 );
//     for( int i = 0; i < number; i++ ) pointer[i] = value;
// }



// inline void sort_integers( int* start, int length )
// {
//     Assert( start != nullptr && length >= 0 );
//     for( int i = 1; i < length; i++ )
//     for( int j = 1; j < length; j++ )
//         if( start[j-1] > start[j] ) 
//             std::swap( start[j-1], start[j] ); // TODO: requires <utility>
// }









#endif
