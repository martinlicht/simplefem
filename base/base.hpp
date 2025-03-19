#ifndef INCLUDEGUARD_BASIC_HPP
#define INCLUDEGUARD_BASIC_HPP


#if __cplusplus < 201402L
#error Compilation of this software requires at least C++14.
#endif



// ============================================================================
// Minimal standard header inclusions
// ============================================================================

#include <cmath>        // pow,
#include <cstdint>      // std::unint_max

#include <limits>       // for numeric_limits
#include <string>       
#include <type_traits>  // for template requirements


// ============================================================================
// Include header file for debug macros
// ============================================================================

#include "debug.hpp"


// ============================================================================
// Define important macros
// ============================================================================

#ifdef ELIDE_HOT_FUNCTIONS
#if defined(__GNUC__) or defined(__clang__)
#define HOTCALL __attribute__((hot,warning("Performance-critical function call not elided.")))
#endif // defined(__GCNUC__) or defined(__clang__)
#else 
#define HOTCALL 
// #warning No hot calls
#endif // ELIDE_HOT_FUNCTIONS


#if defined(__GNUC__) or defined(__clang__)
#define PACKED 
// #define PACKED __attribute__((packed)) // DEACTIVATED BECAUSE IT PREVENTS REFERENCING IN CLANG
#else 
#define PACKED
#endif 


#if __cplusplus >= 202002L
#define LIKELY   [[likely]]
#define UNLIKELY [[unlikely]]
#else
#define LIKELY
#define UNLIKELY
#endif


#if defined(__GNUC__) || defined(__clang__)
#define UNUSED __attribute__((unused))
#else
#define UNUSED
#endif








// ============================================================================
// Casting-type for intentionally discarded variables
// ============================================================================

typedef void void_discard;




// ============================================================================
// Floating-point definitions
// ============================================================================

#if defined(EXTENDED_PRECISION) && defined(SINGLE_PRECISION)
#error Cannot request extended and single precision at the same time!
#endif

#if defined(EXTENDED_PRECISION)
typedef long double Float;
#elif defined(SINGLE_PRECISION)
typedef float Float;
#else 
typedef double Float;
#endif


template<typename T>
constexpr typename std::enable_if< std::is_floating_point<T>::value, T>::type Sqrt( T a, int i = 40 )
{
    T x = a;
    if( not ( x > 0. ) ) return 0.;
    while ( i --> 0 ) x = ( x + a / x ) / 2.f;
    return x;    
}


static const constexpr Float notanumber = std::numeric_limits<Float>::quiet_NaN();

static const constexpr Float machine_epsilon = std::numeric_limits<Float>::epsilon();

static const constexpr Float desired_precision = 
                                    sizeof(Float) == sizeof(float) ? 1e-5 : Sqrt( machine_epsilon );

static const constexpr Float desired_closeness = 
                                    sizeof(Float) == sizeof(float) ? 1e-5 : Sqrt( machine_epsilon );

static const constexpr Float desired_closeness_for_sqrt = 
                                    sizeof(Float) == sizeof(float) ? 1e-5 : 100 * desired_closeness;








// ============================================================================
// Char and string constants
// ============================================================================

static const constexpr char space = ' ';

static const constexpr char* emptystring = "";

static const constexpr char nl = '\n';

static const constexpr char tab = '\t';














// ============================================================================
// System setup
// 
// We use a global singleton for any relevant setup and initial outputs,
// which is called before the main function
// ============================================================================

struct SystemSetup final
{
    SystemSetup() noexcept;
    ~SystemSetup() noexcept;
};

extern const SystemSetup system_setup;









// ============================================================================
// use this to safely cast size_types to C++ int
// ============================================================================

inline constexpr int SIZECAST( std::uintmax_t size )
{
    Assert( size < std::numeric_limits<int>::max() );
    return static_cast<int>( size );
}

// ============================================================================
// use this to safely get C array lengths
// ============================================================================

template < class T, size_t N >
constexpr size_t countof( const T (&array)[N] ) {
  return N;
}




// ============================================================================
// Simple auxiliary arithmetics
// ============================================================================

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
    if( std::is_integral<T>::value && std::is_signed<T>::value )
        assert( x != std::numeric_limits<T>::min() );
    if( x >= 0 ) return  x;
    if( x <= 0 ) return -x;
    Assert( not std::isfinite(x) );
    return x;
}

template<typename T>
inline constexpr T sign( const T& x )
{
    // Assert( std::isfinite(x) );
    if( x > 0 ) return  1;
    if( x < 0 ) return -1;
    else        return  0;
}

template<typename T>
inline constexpr int sign_integer( const T& x )
{
    // Assert( std::isfinite(x) );
    if( x > 0 ) return  1;
    if( x < 0 ) return -1;
    else        return  0;
}

/*
template<typename T>
inline constexpr T maximum( const T& a, const T& b )
{
    if( a >= b ) return a;
    if( a <= b ) return b;
    Assert( ( not std::isfinite(a) ) or ( not std::isfinite(b) ) );
    return a;
//     Assert( a >= b or a <= b ); if( a >= b ) return a; else return b;
}
*/

template<typename T>
inline constexpr T maximum( const T& a )
{
    return a;
}

template<typename T, typename... Args >
inline constexpr T maximum( T a, Args... args )
{
    const T& b = maximum( args... );
    if( a >= b ) return a;
    if( a <= b ) return b;
    Assert( ( not std::isfinite(a) ) or ( not std::isfinite(b) ) );
    return a;
//     Assert( a >= b or a <= b ); if( a >= b ) return a; else return b;
}


template<typename T>
inline constexpr T maxabs( const T& a )
{
    return absolute(a);
}

template<typename T, typename... Args >
inline constexpr T maxabs( T a, Args... args )
{
    const T& b = maxabs( args... );
    return maximum( absolute(a), absolute(b) );
}


template<typename T>
inline constexpr T minimum( const T& a )
{
    return a;
}

template<typename T, typename... Args >
inline constexpr T minimum( T a, Args... args )
{
    const T& b = minimum( args... );
    if( a >= b ) return b;
    if( a <= b ) return a;
    Assert( ( not std::isfinite(a) ) or ( not std::isfinite(b) ) );
    return b;
//     Assert( a >= b or a <= b ); if( a >= b ) return a; else return b;
}

/*
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
    return maximum( t, static_cast<T>( maximum( args... ) ) );
}

template<typename T, typename... Args >
inline constexpr T maxabs( T t, Args... args )
{
    return maxabs( t, static_cast<T>( maxabs( args... ) ) );
}

template<typename T, typename... Args >
inline constexpr T minimum( T t, Args... args )
{
    return minimum( t, static_cast<T>( minimum( args... ) ) );
}
*/


template<typename T>
inline constexpr T square( const T& x )
{
    return x * x;
}

inline constexpr bool is_numerically_small( Float value, Float threshold = desired_closeness )
{
    return absolute(value) < threshold;
}

inline constexpr bool is_numerically_small_sqrt( Float value, Float threshold = desired_closeness_for_sqrt )
{
    return absolute(value) < threshold;
}

inline constexpr bool is_numerically_close( Float value1, Float value2, Float threshold = desired_closeness )
{    
    if( std::isinf(value1) or std::isinf(value2) ) return value1 == value2;
    assert( std::isfinite(value1) and std::isfinite(value2) );
    return is_numerically_small( value1 - value2, threshold );
}

inline constexpr bool is_numerically_one( Float value, Float threshold = desired_closeness )
{
    return is_numerically_close( value, 1., threshold );
}

// https://codingnest.com/the-little-things-comparing-floating-point-numbers/






// ============================================================================
// Power functions
// ============================================================================

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

inline constexpr int power_of_two( int exponent )
{
    assert( exponent >= 0 );
    return 1 << exponent; // return 2 * power_integer( 2, exponent );
}

inline constexpr int sign_power( int exponent )
{
    assert( exponent >= 0 );
    return exponent % 2 == 0 ? 1. : -1;
}





// ============================================================================
//
// Integral factorial, binomials, and auxiliaries
//
// NOTE:
// For small inputs, the naive method seems to perform best,
// the table method is only slightly slower, 
// and the loop method is consistently slowest.
// The differences are in the range of 5%, so fairly small for practical purposes.
//   
// ============================================================================

/*
 * Recursively divide the integer n by larger and larger numbers 1, 2, 3, ... without remainder
 * until the divisor is larger than n. That divisor is the largest numbers whose factorial is at most n.
 */

template<typename T>
inline constexpr uintmax_t largest_factorial_base_AUX( T n, intmax_t k )
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





inline constexpr uintmax_t factorial_integer_table_old( intmax_t n )
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
        default: impossible();
    }
    impossible();
}

inline constexpr uintmax_t factorial_integer_table( intmax_t n )
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

inline constexpr uintmax_t factorial_integer_naive( intmax_t n )
{
    Assert( 0 <= n and n <= 20 );
    if( n == 0 ) { 
        return 1;
    } else {
        return n * factorial_integer_naive( n-1 );
    }
}

inline constexpr uintmax_t factorial_integer_loop( intmax_t n )
{
    Assert( 0 <= n and n <= 20 );
    uintmax_t ret = 1;
    while( n > 0 ) ret *= n--;
    return ret;
}






inline constexpr int factorial_integer( intmax_t n )
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

inline constexpr int binomial_integer( intmax_t n, intmax_t k )
{
    Assert( 0 <= n, "Negative n for integer binomial: ", n ); 
    Assert( 0 <= n );
    
    if( k < 0 or n < k ) return 0;

    uintmax_t result = factorial_integer(n) / ( factorial_integer(k) * factorial_integer(n-k) );
    
    Assert( result <= std::numeric_limits<int>::max() );
    return static_cast<int>(result);
}

inline unsigned long long binomial_integer_secured( unsigned int n, unsigned int k )
{
    Assert( k <= n ); 

    if( k > n - k ) { 
        k = n - k;
    }
    
    unsigned long long result = 1;

    for( unsigned int i = 1; i <= k; i++ ) 
    {
        // Check for overflow before multiplication
        Assert( result <= std::numeric_limits<unsigned long long>::max() / (n - i + 1) );
        
        // Multiply result by (n - i + 1) and then divide result by i (guaranteed to be exact)
        result *= (n - i + 1);
        result /= i;
    }

    return result;
}






// ============================================================================
//
// Numerical factorial, binomials, and auxiliaries
//
// NOTE:
// For small inputs, the loop method seems to be fastest, 
// whereas the recursive form of the factorial performs 10% slower.
//   
// ============================================================================

inline constexpr Float factorial_numerical_naive( intmax_t n )
{
    Assert( 0 <= n );
    if( n == 0 ) { 
        return 1.;
    } else {
        return static_cast<Float>(n) * factorial_numerical_naive( n-1 );
    }
}

inline constexpr Float factorial_numerical_loop( intmax_t n )
{
    Assert( 0 <= n );
    Float ret = 1.;
    while( n > 0 ) ret *= static_cast<Float>(n--);
    return ret;
}

inline constexpr Float factorial_numerical_table( intmax_t n )
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

inline constexpr Float factorial_numerical( intmax_t n )
{
    return factorial_numerical_naive( n );
}


inline constexpr Float binomial_numerical( intmax_t n, intmax_t k )
{
    Assert( 0 <= n );
    if( k < 0 or n < k )
        return 0.;
    return factorial_numerical(n) / ( factorial_numerical(k) * factorial_numerical(n-k) );
}



























// ============================================================================
// Time utilities
// ============================================================================

typedef uintmax_t timestamp;

timestamp timestampnow();

std::string timestamp2measurement( const timestamp& t );

std::string measurementnow();

std::string timestamp2digitalcode( const timestamp& t );

std::string digitalcodenow();
















// ============================================================================
// String utilities
// ============================================================================

// ----------------------------------------------------------------------------
// count the white space within STL string
// ----------------------------------------------------------------------------

int count_white_space( const std::string& str ); 

// ----------------------------------------------------------------------------
// insert tabs before each line
// ----------------------------------------------------------------------------

std::string tab_each_line( std::string str );

// ----------------------------------------------------------------------------
// convert a C string to an integer
// ----------------------------------------------------------------------------

int string_to_integer( const char* s, const char* __restrict__ *endptr, unsigned int base, bool& has_overflown );

// TODO(martinlicht): These functions are not used in the base. Move to utilities 














// ============================================================================
// Generic streaming
// ============================================================================

// ----------------------------------------------------------------------------
// printf into C++ strings and streams 
// ----------------------------------------------------------------------------

std::string printf_into_string( const char* formatstring, ... ) 
__attribute__ (( format (printf,1,2) ));

#if __cplusplus < 202002L
#define printf_into_stream( stream, ... ) { stream << printf_into_string( __VA_ARGS__ ); }
#else
#define printf_into_stream( stream, formatstring, ... ) { stream << printf_into_string( formatstring __VA_OPT__(,) __VA_ARGS__ ); }
#endif
// template< typename L, typename... Params >
// inline void printf_into_stream( L& stream, const char* formatstring, Params... args )
// {
//     stream << printf_into_string( formatstring, args... );
// }










#endif
