
#include <cctype>   // char is ... checks
#include <cfenv>    // For controlling the floating-point behavior
#include <cfloat>   // FLT_EVAL_METHOD, _controlfp_s on Windows 
#include <cstdio>   // print
#include <cstdarg>  // variadic function macros
#include <ctime>    // localtime 

#include <chrono>
#include <limits> // for std::float_round_style
#include <string>
#include <type_traits>

#include "base.hpp"
#include "logging.hpp"


#if defined(__SSE__)
#include <xmmintrin.h> // _MM_SET_FLUSH_ZERO_MODE
#endif
#if defined(__SSE3__)
#include <pmmintrin.h> // _MM_SET_DENORMALS_ZERO_MODE
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif // #if defined(_OPENMP)


// ===============================================================================================
// The following global variable (!) will be constructed before the control flow enters 
// the main function. Its constructor will print several system settings on the standard output
// and additionally set several floating-point switches.
// ===============================================================================================

const SystemSetup system_setup;

SystemSetup::SystemSetup() noexcept
{
    // ===============================================================================================
    // Print numerous information on the compilation configuration:
    // - Current git commit 
    // - _DEBUG (MSVC, whether debug is on)
    // - NDEBUG
    // - __OPTIMIZE__ (GCC< whether any optimization is on)
    // - Operating system
    // ===============================================================================================
    
    #if defined(GIT_COMMIT_ID)
    LOG << "---\tCurrent Git commit ID: " << GIT_COMMIT_ID << nl;
    #endif

    const time_t t = std::time(nullptr);
    const tm zeit = *std::localtime(&t); // DONTFIX: not thread-safe but this function will only be executed once.
    LOGPRINTF("---\t%d-%02d-%02d %02d:%02d:%02d\n", zeit.tm_year + 1900, zeit.tm_mon + 1, zeit.tm_mday, zeit.tm_hour, zeit.tm_min, zeit.tm_sec );
    
    LOGPRINTF("---\tCompiler version: %s\n", __VERSION__);

    #ifdef _DEBUG
    LOGPRINTF("---\tMSVC Debugging flags enabled.\n");
    #else
    LOGPRINTF("---\tMSVC Debugging flags not enabled.\n");
    #endif

    #ifdef NDEBUG
    LOGPRINTF("---\tNDEBUG enabled.\n");
    #else
    LOGPRINTF("---\tNDEBUG not enabled.\n");
    #endif

    #ifdef __OPTIMIZE__
    LOGPRINTF("---\tCompiler optimization level: %d\n", __OPTIMIZE__);
    // Add more performance-related information here
    #endif

    #ifdef __unix__
    LOGPRINTF("---\tOS: Unix\n");
    #elif defined(_WIN32)
    LOGPRINTF("---\tOS: Windows\n");
    #else
    #error "Unknown platform"
    #endif


    // ===============================================================================================
    // Print information on OpenMP, if enabled.
    // ===============================================================================================
    
    #if defined(_OPENMP)
    
    LOG << "---\tOpenMP Value                   : " << _OPENMP                                                                     << nl;
    LOG << "---\tMaximum number of threads      : " << omp_get_max_threads()                                                       << nl;
    LOG << "---\tMaximum active levels          : " << omp_get_max_active_levels()                                                 << nl;
    LOG << "---\tThread limit                   : " << omp_get_thread_limit()                                                      << nl;
    LOG << "---\tNumber of processors available : " << omp_get_num_procs()                                                         << nl;
    LOG << "---\tDefault thread affinity        : " << ( omp_get_proc_bind() == omp_proc_bind_false ? "No affinity" : "Affinity" ) << nl;
    LOG << "---\tmax num nested active parallel : " << omp_get_max_active_levels()                                                 << nl;
    // omp_get_nested is deprecated
    // LOG << "---\tNested parallelism supported:   " << ( omp_get_nested() ? "Yes" : "No" )                                         << nl;
    LOG << "---\tMaximum number of places       : " << omp_get_num_places()                                                        << nl;
    LOG << "---\tDynamic adjustment of the number of threads enabled: " << ( omp_get_dynamic() ? "Yes" : "No" ) << nl;
    for( int p = 0; p < omp_get_num_places(); p++ ) {
        LOG << "---\t\tMaximum number of processors per place: " << p << " - > " << omp_get_place_num_procs(p) << nl;
    }
    
    // // provoke the spawning of the OpenMP threads with some parallel block
    // #pragma omp parallel
    // {
    //     printf("<%d/%d>", omp_get_thread_num(), omp_get_num_threads() );
    // }
    // printf("\n");
    
    #else
    LOG << "---\tOpenMP is not enabled.\n";
    #endif

    // #ifdef _OPENMP
    // LOGPRINTF("OpenMP is enabled.\n");
    // LOGPRINTF("Number of threads: %d\n", omp_get_max_threads());
    // LOGPRINTF("Nested parallelism supported: %s\n", omp_get_nested() ? "Yes" : "No");
    // LOGPRINTF("Dynamic adjustment of the number of threads enabled: %s\n", omp_get_dynamic() ? "Yes" : "No");
    // LOGPRINTF("Number of processors available: %d\n", omp_get_num_procs());
    // LOGPRINTF("Default thread affinity: %s\n", omp_get_proc_bind() == omp_proc_bind_false ? "No affinity" : "Affinity");
    // #else
    // LOGPRINTF("OpenMP is not enabled.\n");
    // #endif
    
    
    // ===============================================================================================
    // Print numerous information on the floating-point environment.
    // ===============================================================================================
    
    LOG << "---\tFLT_EVAL_METHOD: " << FLT_EVAL_METHOD << space;

    switch (FLT_EVAL_METHOD) {
        case -1:
            LOG << "The default precision is not known." << nl;
            break;
        case 0:
            LOG << "All operations and constants evaluate in the range and precision of the type used." << nl;
            break;
        case 1:
            LOG << "All operations and constants evaluate in the range and precision of double." << nl;
            break;
        case 2:
            LOG << "All operations and constants evaluate in the range and precision of long double." << nl;
            break;
        default:
            if (FLT_EVAL_METHOD < 0) {
                LOG << "Implementation-defined behavior." << nl;
            } else {
                LOG << "Unknown evaluation mode." << nl;
            }
            break;
    }


    LOG << "---\tFLT_ROUNDS: " << FLT_ROUNDS << space;

    switch (FLT_ROUNDS) {
        case std::round_indeterminate:
            LOG << "Rounding style cannot be determined." << nl;
            break;
        case std::round_toward_zero:
            LOG << "Rounding toward zero." << nl;
            break;
        case std::round_to_nearest:
            LOG << "Rounding toward nearest representable value." << nl;
            break;
        case std::round_toward_infinity:
            LOG << "Rounding toward positive infinity." << nl;
            break;
        case std::round_toward_neg_infinity:
            LOG << "Rounding toward negative infinity." << nl;
            break;
        default:
            LOG << "Unknown rounding style." << nl;
            break;
    }


    #if __cplusplus >= 201703L
    
    LOG << "---\tFLT_HAS_SUBNORM: " << FLT_HAS_SUBNORM << space;
    switch (FLT_HAS_SUBNORM) {
        case -1:
            LOG << "Support for subnormal numbers in float is indeterminable." << nl;
            break;
        case 0:
            LOG << "Float type does not support subnormal numbers." << nl;
            break;
        case 1:
            LOG << "Float type supports subnormal numbers." << nl;
            break;
        default:
            LOG << "Unknown value for FLT_HAS_SUBNORM." << nl;
            break;
    }

    LOG << "---\tDBL_HAS_SUBNORM: " << DBL_HAS_SUBNORM << space;
    switch (DBL_HAS_SUBNORM) {
        case -1:
            LOG << "Support for subnormal numbers in double is indeterminable." << nl;
            break;
        case 0:
            LOG << "Double type does not support subnormal numbers." << nl;
            break;
        case 1:
            LOG << "Double type supports subnormal numbers." << nl;
            break;
        default:
            LOG << "Unknown value for DBL_HAS_SUBNORM." << nl;
            break;
    }

    LOG << "---\tLDBL_HAS_SUBNORM: " << LDBL_HAS_SUBNORM << space;
    switch (LDBL_HAS_SUBNORM) {
        case -1:
            LOG << "Support for subnormal numbers in long double is indeterminable." << nl;
            break;
        case 0:
            LOG << "Long double type does not support subnormal numbers." << nl;
            break;
        case 1:
            LOG << "Long double type supports subnormal numbers." << nl;
            break;
        default:
            LOG << "Unknown value for LDBL_HAS_SUBNORM." << nl;
            break;
    }
    
    #endif 

    #if defined(_WIN32)
    {
        // Query the current floating-point control word
        unsigned int currentControlWord = _controlfp(0, 0);
        
        // Check precision control bits
        switch(currentControlWord & _MCW_PC) {
            case _PC_24: 
                LOG << "---\tPrecision control: 24-bit." << nl;
                break;
            case _PC_53:
                LOG << "---\tPrecision control: 53-bit." << nl;
                break;
            case _PC_64:
                LOG << "---\tPrecision control: 64-bit." << nl;
                break;
            default:
                LOG << "---\tPrecision control: unknown." << nl;
                break;
        }
    }    
    #endif


    // ===============================================================================================
    // The following settings enable the flushing of denormals, both operands and results, to zero.
    // 
    // If OpenMP is active, then we execute it in a parallel region, 
    // so that each thread sets the flushing mode 
    //  1. We assume that the main thread participates in the parallel regions, and that
    //  2. the threads are persistent between parallel regions
    // If OpenMP is active, we notify on the main thread
    // 
    // Generally, we use the floating-point intrinsics headers. 
    // On Windows, however, we use the native headers
    // ===============================================================================================
    
    #if defined(_OPENMP)
    #pragma omp parallel
    {
        if( omp_get_thread_num() == 0 )
    #endif    

    #if defined(_WIN32)
        LOGPRINTF("---\tFlushing subnormal numbers\n");
        _controlfp_s( nullptr, _DN_FLUSH, _MCW_DN );        // Use the native Windows interface
    #elif defined(__SSE__)
        LOGPRINTF("---\tFlushing subnormal numbers\n");     // xmmintrin.h
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
        #ifdef __SSE3__
            _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON); // pmmintrin.h
        #endif
    #endif

    #if defined(_OPENMP)
        (void)0; }
    #endif    

}

SystemSetup::~SystemSetup() noexcept
{
    // #if defined(_OPENMP)
    // LOG << "---\tSystem Reporter: finished\n";
    // #endif
}   










// ===============================================================================================
// Text-related functions
// ===============================================================================================


// Since all floating-point literals throughout are double unless marked otherwise 
// we enforce that `Float` is at least enough to store double. Any of those should do:
// 
// static_assert( Float(std::numeric_limits<double>::max()) == std::numeric_limits<double>::max(), "Float must be at least double" );
// static_assert( sizeof(Float) >= sizeof(double), "Float must be at least double" );

int count_white_space( const std::string& str ) 
{ 
    int ret = 0;  
    for( int c = 0; c < str.size(); c++ ) if( 0 != std::isspace( str[c] ) ) ret++; 
    return ret;
} 

std::string tab_each_line( std::string str ) 
{ 
    str.insert( 0, 1, '\t' );
    for( int c = str.size(); c > 0; c-- ) {
        if( str[c-1] == '\n' )
            str.insert(c, 1, '\t' );
    }
    return str;
} 


int string_to_integer( const char* s, const char* __restrict__ *endptr, unsigned int base, bool& has_overflown ) {
    
    assert( s != nullptr );
    assert( ( 2 <= base and base <= 36 ) or ( base == 0 ) );

    unsigned int result = 0;
    
    int sign = 1;

    has_overflown = false;
    
    // Skip leading whitespace
    while( std::isspace( (unsigned char)*s ) != 0 ) { s++; }

    if( *s == '\0' ) {
        // printf("Only white space\n");
        return 0;
    } 
        
    // Handle optional sign
    if( *s == '-' ) { 
        sign = -1;
        s++;
    } else if( *s == '+' ) {
        s++;
    }

    if( *s == '\0' ) { 
        // printf("Only signs\n");
        return 0;
    }
        
    // Auto-detect if base equals zero 
    if( base == 0 ) {
        
        // printf("Auto-detect base\n");
        
        base = 10;

        if( *s == '0' ) { 
            
            base = 8;
            s++;
        
            if( *s == '\0' ) return 0;
            
            if( *s == 'x' || *s == 'X' ) {
                base = 16; 
                s++;
            }

        }
    }

    // Convert digits
    while( *s != '\0' ) {
        
        int digit;
        
        // Find digits within the correct range 
        if( *s >= '0' && *s <= '9' ) {
            digit = *s - '0';
        } else if( *s >= 'a' && *s <= 'z' ) {
            digit = *s - 'a' + 10;
        } else if( *s >= 'A' && *s <= 'Z' ) {
            digit = *s - 'A' + 10;
        } else {
            break;
        }

        // printf("Digit detected: %u\n", digit );

        if( digit >= base ) break;

        if( not has_overflown and result > ( std::numeric_limits<int>::max() - digit) / base) {
            // printf("Overflowing... %u %u\n", result, digit );
            has_overflown = true;
        } else {
            result = result * base + (unsigned int)digit;
        }

        s++;
    }

    // Store pointer one past the last digit 
    if( endptr != nullptr ) 
    {
        *endptr = s;
    }

    // Assuming two's complement (C++20 and C23)
    if( sign == -1 and result > 1u + (unsigned int)(std::numeric_limits<int>::max()) ) has_overflown = true;
    
    if( has_overflown ) 
    {
        // printf("Overflown\n");
        return ( sign == 1 ) ? ( std::numeric_limits<int>::max() ) : ( std::numeric_limits<int>::min() );
    }

    return sign * result;
}














































std::string printf_into_string( const char* formatstring, ... )
{
    
    va_list args;
    
    va_start( args, formatstring );
    const std::size_t length = std::vsnprintf(nullptr, 0, formatstring, args ) + 1;
    va_end( args );
    
    char* c_str = new char[length];
    
    va_start( args, formatstring );
    (void_discard)std::vsnprintf( c_str, length, formatstring, args );
    va_end( args );
    
    std::string ret( c_str );
    delete[] c_str;
    return ret;
}

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





static_assert( std::is_integral< decltype( std::chrono::time_point_cast< std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count() ) >::value , "Time measurement must be integral" );

timestamp timestampnow()
{
    
    static const timestamp start_timestamp = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    const timestamp                    now = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    Assert( now >= start_timestamp );
    
    return now - start_timestamp;
}



std::string timestamp2measurement( const timestamp& t )
{
    return std::to_string( static_cast<uintmax_t>(t) ) + "ms";
}

std::string measurementnow()
{
    return timestamp2measurement( timestampnow() );
}



std::string timestamp2digitalcode( const timestamp& t )
{
    const int numdigits = 10;
    const int fulllength = numdigits+1;
    char digits[fulllength];
    assert( t < 9999999999 ); // ca. 115 days 
    (void_discard)snprintf( digits, fulllength, "%*ju", numdigits, (uintmax_t)t );
    for( int i = 0; i < numdigits; i++ ) if( digits[i] == ' ' ) digits[i] = '_';
    return std::string(digits);
}

std::string digitalcodenow()
{
    return timestamp2digitalcode( timestampnow() );
}



// std::string protocolprefixnow()
// {
//     // static const std::string foo = std::string("\e[36m[");
//     // static const std::string bar = std::string("]\e[39m\t");
//     static const std::string foo = std::string("[");
//     static const std::string bar = std::string("]\t");
//     return foo + digitalcodenow() + bar;
// }










