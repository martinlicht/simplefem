#ifndef INCLUDEGUARD_DEBUG_HPP
#define INCLUDEGUARD_DEBUG_HPP


/* Definitions for assert macros 
 * 
 * The general structure of this framework is as follows.
 * We define the macros 
 * 
 *  - unreachable()
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
 */







/* 
 * 
 * Unless we are told to use the original assert macros, we define the internal functions 
 * that implement the error output. No checks are performed here!
 * 
 * If the original assert macro is to be used, then we don't need those.
 * 
 */

#ifndef USE_ORIGINAL_ASSERT_MACRO

#include <cstdio>
#include <cstdlib>

inline void myActualAssert [[noreturn]] ( const char* filename, const int linenumber, const char* expression, const char* message )
{
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );
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
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );    
#ifdef __cpp_exceptions
    throw(0);
#else // __cpp_exceptions
    abort();
#endif // __cpp_exceptions    
}

inline void myActualUnreachable [[noreturn]] ( const char* filename, const int linenumber )
{
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "!!\n" ); 
    fprintf( stderr, "!!\tUnreachable code reached:\n!!!!\t%s:%d\n", filename, linenumber ); 
    fprintf( stderr, "!!\n" ); 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
#ifdef __cpp_exceptions
    throw(0);
#else // __cpp_exceptions
    abort();
#endif // __cpp_exceptions    
}

inline void myActualUnimplemented [[noreturn]] ( const char* filename, const int linenumber )
{
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "!!\n" ); 
    fprintf( stderr, "!!\tUnimplemented execution path reached:\n!!!!\t%s:%d\n", filename, linenumber ); 
    fprintf( stderr, "!!\n" ); 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ); 
#ifdef __cpp_exceptions
    throw(0);
#else // __cpp_exceptions
    abort();
#endif // __cpp_exceptions    
}

#endif // USE_ORIGINAL_ASSERT_MACRO








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

// The following contains the framework to enable a varyadic assert macro
// The internal function is templated; after the first few standard arguments
// all remaining arguments are put into a templated function 
// that concatenates those arguments into a string. 
// If 
//   - no extra arguments are there, an empty string is produced 
//   - there are extra arguments, they concatenated into a string, with separators
// 
// The stringification uses the shift operator into a stringstream


#include <string>
#include <sstream>

// nothing to concat: empty string
inline std::string Concat2String()
{
    return "";
}

// one argument to stringify, base of induction 
template< typename T >
inline std::string Concat2String( const T& t )
{
    std::ostringstream ss;
    ss << t;
    return ss.str();
}

// recursively build a string by stringifying arguments,
// with separators in between. Base case has only one argument.
template< typename T, typename... Params >
inline std::string Concat2String( const T& t, const Params&... params )
{
    std::ostringstream ss;
    ss << t << '\t' << Concat2String( params... );
    return ss.str();
}

#endif // !defined NDEBUG && !defined DISCARD_ASSERT_MESSAGES && !defined USE_ORIGINAL_ASSERT_MACRO 



/* 
 * 
 * Now we define the relevant macros. 
 * If assertions are disabled, then `Assert` does nothing 
 * and `unreachable` and `unimplemented` simply call exit with failure.
 * 
 * 
 * Now, if we use the original assert macro, then `Assert` simply redirects the boolean variable into the standard `assert` macro.
 * If the C++ version is below 20, we wrap this into another varyadic macro, which clogs the output a bit
 * If the C++ version at least 20, then we forward the boolean via a varyadic macro
 * 
 * Instead, if the custom assert macro is used, then 
 * If the C++ version is below 20, we do as above 
 * If the C++ version is least 20, we use everything
 * 
 */



#ifdef NDEBUG

#include <cassert>
#define Assert(...)   (static_cast<void>(0))
#define unreachable()   (static_cast<void>(0),exit( EXIT_FAILURE ))
#define unimplemented() (static_cast<void>(0),exit( EXIT_FAILURE ))

#else // NDEBUG

#ifdef USE_ORIGINAL_ASSERT_MACRO

#include <cassert>

#if __cplusplus < 202002L
constexpr bool Cond( bool b ) { return b; }
template<typename... Params> constexpr bool Cond( bool b, const Params&... params ) { return b; }
#define Assert(...)     assert(Cond(__VA_ARGS__))
#else
#define Assert(x,...)     assert(x)
#endif // __cplusplus < 202002L

#define unreachable()   fprintf( stderr, "Unreachable reached: %s:%d\n", __FILE__, __LINE__ ),abort()
#define unimplemented() fprintf( stderr,  "Unimplemented path: %s:%d\n", __FILE__, __LINE__ ),abort()

#else // USE_ORIGINAL_ASSERT_MACRO

#if __cplusplus < 202002L
#include <cassert>
constexpr bool Cond( bool b ) { return b; }
template<typename... Params> constexpr bool Cond( bool b, const Params&... params ) { return b; }
#define Assert(...)     assert(Cond(__VA_ARGS__))
#else 
#ifndef DISCARD_ASSERT_MESSAGES
#define Assert(x,...) (static_cast<bool>(x)?(void(0)):myActualAssert( __FILE__, __LINE__, #x, 0 __VA_OPT__(+1)?Concat2String(__VA_ARGS__).c_str():nullptr ) )
#else
#define Assert(x,...) (static_cast<bool>(x)?(void(0)):myActualAssert( __FILE__, __LINE__, #x, nullptr ) )
#endif // DISCARD_ASSERT_MESSAGES
#endif // __cplusplus < 202002L


#define unreachable()   { myActualUnreachable(__FILE__, __LINE__), abort(); }
#define unimplemented() { myActualUnimplemented(__FILE__, __LINE__), abort(); }

#endif //USE_ORIGINAL_ASSERT_MACRO

#endif //NDEBUG


#endif //INCLUDEGUARD_DEBUG_HPP
