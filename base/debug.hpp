#ifndef INCLUDEGUARD_BASE_DEBUG_HPP
#define INCLUDEGUARD_BASE_DEBUG_HPP


/* Definitions for assert macros 
 * 
 * The general structure of this framework is as follows.
 * We define the macros 
 * 
 *  - impossible()
 *  - unimplemented()
 *  - Assert(x)
 *  - Assert(x,...)
 * 
 * Those definitions are filled up in different ways.
 * 
 * 1)
 * If USE_ORIGINAL_ASSERT_MACRO is set, then we use 
 * the capabilities of the C library, in particular the
 * original `assert` macro.
 * 
 * Otherwise, we define them by ourselves.
 * 
 * 2)
 * There is another case distinction depending on whether 
 * we have set NDEBUG or not. If NDEBUG is set, then we 
 * fix the terms as empty. Otherwise non-trivial definitions
 * follow.
 * 
 * 
 * The macro `unreachable` clashes with the C++23 unreachable
 * 
 * 
 */




// ============================================================================
// Define a few internal macros
// ============================================================================

/* 
 * For the unreachable/impossible macro, we may use a built-in signal to mark the code unreachable,
 * provided that this is available on this platform 
 */

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define UNREACHABLE_SIGNAL() __builtin_unreachable()
#elif defined(_MSC_VER)
#define UNREACHABLE_SIGNAL() __assume(0)
#else
#define UNREACHABLE_SIGNAL() 
#endif

/* 
 * Irrespective of whether we use it or not, we define a termination command 
 */
#ifdef __cpp_exceptions
#define TERMINATION_COMMAND() throw(0)
#else // __cpp_exceptions
#define TERMINATION_COMMAND() abort()
#endif // __cpp_exceptions    




// ============================================================================
// Unless the original assert macro is defined,
// we define the internal functions for debugging.
// Otherwise, we can leave this out. 
// ============================================================================

/* 
 * 
 * Unless we are told to use the original assert macros, we define the internal functions 
 * that implement the error output. No checks are performed here!
 * 
 * If the original assert macro is to be used, then we don't need those.
 * 
 */

#ifndef USE_ORIGINAL_ASSERT_MACRO

#include <cstdio>   // printf
#include <cstdlib>  // possibly for abort()

// If a backtracer is requested, then provide if available
#if defined(USE_BACKTRACER) 
#include <version>
#if defined(__cpp_lib_stacktrace)
#include <stacktrace> // stack trace capability (C++23)
inline void myStackTrace( unsigned int skip )
{
    unsigned int i_max = 20;
    auto st = std::stacktrace::current();
    fprintf( stderr, "!! printing stack trace (%u-%u):\n", skip, i_max );
    for( std::size_t i = skip; i < 20 && i < st.size(); i++ ) 
    {
        const auto& entry = st.at(i);
        fprintf( stderr, "!! #%5lu %s\n",   i, entry.description().c_str() );
        fprintf( stderr, "!!      @ %s:%u \n", entry.source_file().c_str(), entry.source_line() );
    }
    fprintf( stderr, "-------------------------------\n");
}
#else
#error "Backtracer requested but unavailable."
#endif
#endif
// #endif


static void myActualUnreachable [[noreturn]] ( const char* filename, const int linenumber )
{
    fprintf( stderr, "\n" ); 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "!!\n" ); 
    fprintf( stderr, "!!\tUnreachable code reached:\n" ); 
    fprintf( stderr, "!!\t%s:%d\n", filename, linenumber ); 
    fprintf( stderr, "!!\n" ); 
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );

    #if defined(USE_BACKTRACER) && defined(__cpp_lib_stacktrace)
    myStackTrace(2); // Print a stack trace if requested and available
    #endif

    TERMINATION_COMMAND();
}

inline void myActualUnimplemented [[noreturn]] ( const char* filename, const int linenumber )
{
    fprintf( stderr, "\n" ); 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "!!\n" ); 
    fprintf( stderr, "!!\tUnimplemented execution path reached:\n" ); 
    fprintf( stderr, "!!\t%s:%d\n", filename, linenumber ); 
    fprintf( stderr, "!!\n" ); 
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 

    #if defined(USE_BACKTRACER) && defined(__cpp_lib_stacktrace)
    myStackTrace(2); // Print a stack trace if requested and available
    #endif

    TERMINATION_COMMAND();
}

inline void myActualAssert [[noreturn]] ( const char* filename, const int linenumber, const char* expression, const char* message )
{
    fprintf( stderr, "\n" ); 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "!!\n" );
    fprintf( stderr, "!!\tThe following assertion failed.\n" );
    fprintf( stderr, "!!\t%s,l.%d\n", filename, linenumber );
    fprintf( stderr, "!!\n" );
    fprintf( stderr, "!!\t\t%s\n", expression );
    if( message != nullptr ){
    fprintf( stderr, "!!\n" );
    fprintf( stderr, "!!\t\t%s\n", message );
    }
    fprintf( stderr, "!!\n" );
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );

    #if defined(USE_BACKTRACER) && defined(__cpp_lib_stacktrace)
    myStackTrace(2); // Print a stack trace if requested and available
    #endif

    TERMINATION_COMMAND();
}

#endif // USE_ORIGINAL_ASSERT_MACRO







// ============================================================================
// Compile-time string concatenation (if necessary for debug build)
// ============================================================================

/* 
 * 
 * Unless 
 * 
 * - assertions are disabled 
 * OR
 * - we discard assert messages 
 * OR 
 * - we use the original assert macro
 * 
 * the following prepares facilities for a varyadic assert macro: 
 * we provide a varyadic template (C++11) which builds a string from a number of objects.
 * The string is via the shift operator << and assumes that this can be performed. 
 * No checks done here!.
 * 
 * 
 */

#if !defined NDEBUG && !defined DISCARD_ASSERT_MESSAGES && !defined USE_ORIGINAL_ASSERT_MACRO 

// We introduce compile-time string concation that will enable a varyadic assert macro.
// The internal function is templated; after the first few standard arguments.
// all remaining arguments are put into a templated function 
// that concatenates those arguments into a string. 
// 
//   - If no extra arguments are there, an empty string is produced 
//   - If, instead, there are extra arguments, they concatenated into a string, with separators
// 
// The stringification uses the shift operator into a stringstream.

#include <string>
#include <sstream>

// nothing to concat: empty string
inline std::string Concat2String()
{
    return "";
}

// // one argument to stringify, base of induction 
// template< typename T >
// inline std::string Concat2String( const T& t )
// {
//     std::ostringstream ss; ss.flags( std::ios::scientific ); ss.precision( 15 );
//     ss << t;
//     return ss.str();
// }

// recursively build a string by stringifying arguments,
// with separators in between. Base case has only one argument.
template< typename T, typename... Params >
inline std::string Concat2String( const T& t, const Params&... params )
{
    std::ostringstream ss; ss.flags( std::ios::scientific ); ss.precision( 15 );
    ss << t << '\t' << Concat2String( params... );
    return ss.str();
}

#endif // !defined NDEBUG && !defined DISCARD_ASSERT_MESSAGES && !defined USE_ORIGINAL_ASSERT_MACRO 





// ============================================================================
// Main definitions for users 
// ============================================================================

/* 
 * 
 * Now we define the relevant macros. 
 * If assertions are disabled, then `Assert` does nothing 
 * and `unreachable` and `unimplemented` simply call exit with failure.
 * 
 * 
 * Now, if we use the original assert macro, then `Assert` simply redirects the boolean variable into the standard `assert` macro.
 * If the C++ version is below 20, we wrap this into another varyadic macro, which pollutes the output somewhat
 * If the C++ version at least 20, then we forward the boolean via a varyadic macro
 * 
 * Instead, if the custom assert macro is used, then 
 * If the C++ version is below 20, we do as above 
 * If the C++ version is least 20, we use everything
 * 
 */

#ifdef NDEBUG

/* If debugging is disabled, then the macros are no-ops */

#include <cassert>
#define impossible()    { UNREACHABLE_SIGNAL();  }
#define unimplemented() { UNREACHABLE_SIGNAL();  }
#define Assert(...)     { (static_cast<void>(0));}

// Note: replaced `exit( EXIT_FAILURE ))` with `abort()`. The latter function is best for irregular termination.

#else // NDEBUG

/* We can assume that debugging is enabled */

#ifdef USE_ORIGINAL_ASSERT_MACRO

/* Now assume the case where we default to the original assert macro */

#include <cassert>

#define impossible()    { fprintf( stderr, "Unreachable reached: %s:%d\n", __FILE__, __LINE__ ); UNREACHABLE_SIGNAL(); }
#define unimplemented() { fprintf( stderr,  "Unimplemented path: %s:%d\n", __FILE__, __LINE__ ); UNREACHABLE_SIGNAL(); }

#if __cplusplus < 202002L
constexpr bool Cond( bool b ) { return b; }
template<typename... Params> constexpr bool Cond( bool b, const Params&... params ) { return b; }
#define Assert(...)     assert(Cond(__VA_ARGS__))
#else
#define Assert(x,...)     assert(x)
#endif // __cplusplus < 202002L

#else // USE_ORIGINAL_ASSERT_MACRO

#define impossible()    { myActualUnreachable  (__FILE__, __LINE__); UNREACHABLE_SIGNAL(); }
#define unimplemented() { myActualUnimplemented(__FILE__, __LINE__); UNREACHABLE_SIGNAL(); }

/* Assume the case that debuggig is enabled and we use the custom assert macro           */
/* If the C++ version is before C++20, then we handle assert messages manually           */
/* Else, we handle them depending on whether we are instructed to keep them or toss them */

#if __cplusplus < 202002L
#include <cassert>
inline constexpr bool Cond( bool b ) { return b; }
template<typename T, typename... Params> inline constexpr bool Cond( bool b, const T& param, const Params&... params ) { static_cast<void>(param); return Cond(b,params...); }
#define Assert(...)     assert(Cond(__VA_ARGS__))
#else 
#ifndef DISCARD_ASSERT_MESSAGES
#define Assert(x,...) (static_cast<bool>(x)?(void(0)):myActualAssert( __FILE__, __LINE__, #x, 0 __VA_OPT__(+1)?Concat2String(__VA_ARGS__).c_str():nullptr ) )
#else
#define Assert(x,...) (static_cast<bool>(x)?(void(0)):myActualAssert( __FILE__, __LINE__, #x, nullptr ) )
#endif // DISCARD_ASSERT_MESSAGES
#define assert(x,...) Assert(x,__VA_ARGS__)
#endif // __cplusplus < 202002L


#endif //USE_ORIGINAL_ASSERT_MACRO

#endif //NDEBUG


#endif //INCLUDEGUARD_BASE_DEBUG_HPP
