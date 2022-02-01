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
 * If FLAG_USE_ORIGINAL_ASSERT_MACRO is set, then we use 
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







#include <string>
#include <cstdio>
#include <cstdlib>

inline void myActualAssert( const char* filename, const int linenumber, const char* expression, const std::string message = "" )
{
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );
    fprintf( stderr, "!!\n" );
    fprintf( stderr, "!!\tThe following assertion failed.\n" );
    fprintf( stderr, "!!\t%s,l.%d\n", filename, linenumber );
    fprintf( stderr, "!!\n" );
    fprintf( stderr, "!!\t\t%s\n", expression );
    if( message != "" ){
    fprintf( stderr, "!!\n" );
    fprintf( stderr, "!!\t\t%s\n", message.c_str() );
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
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), 
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), 
    fprintf( stderr, "!!\n" ), 
    fprintf( stderr, "!!\tUnreachable code reached:\n!!!!\t%s:%d\n", __FILE__, __LINE__ ), 
    fprintf( stderr, "!!\n" ), 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), 
#ifdef __cpp_exceptions
    throw(0);
#else // __cpp_exceptions
    abort();
#endif // __cpp_exceptions    
}

inline void myActualUnimplemented [[noreturn]] ( const char* filename, const int linenumber )
{
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), 
    fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), 
    fprintf( stderr, "!!\n" ), 
    fprintf( stderr, "!!\tUnimplemented execution path reached:\n!!!!\t%s:%d\n", __FILE__, __LINE__ ), 
    fprintf( stderr, "!!\n" ), 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), 
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), 
#ifdef __cpp_exceptions
    throw(0);
#else // __cpp_exceptions
    abort();
#endif // __cpp_exceptions    
}









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
    std::stringstream ss;
    ss << t;
    return ss.str();
}

// recursively build a string by stringifying arguments,
// with separators in between. Base case has only one argument.
template< typename T, typename... Params >
inline std::string Concat2String( const T& t, const Params&... params )
{
    std::stringstream ss;
    ss << t << '\t' << Concat2String( params... );
    return ss.str();
}

// template< typename... Params >
// inline void myTemplatedAssert( const char* filename, const int linenumber, const char* expression, const Params&... params )
// {
//     myActualAssert( filename, linenumber, expression, Concat2String( params... ) );
// }








#ifdef NDEBUG

#define Assert(x,...)   (static_cast<void>0)
#define unreachable()   (static_cast<void>0)
#define unimplemented() (static_cast<void>0)

#else // NDEBUG

#ifdef FLAG_USE_ORIGINAL_ASSERT_MACRO

#include <cassert>
#define Assert(x,...)   assert(x)
#define unreachable()   assert(false)
#define unimplemented() assert(false)

#else // FLAG_USE_ORIGINAL_ASSERT_MACRO

#define Assert(x,...) (static_cast<bool>(x)?(void(0)):myActualAssert( __FILE__, __LINE__, #x, Concat2String(__VA_ARGS__) ) )
#define unreachable()   { myActualUnreachable(__FILE__, __LINE__), abort(); }
#define unimplemented() { myActualUnimplemented(__FILE__, __LINE__), abort(); }

#endif //FLAG_USE_ORIGINAL_ASSERT_MACRO

#endif //NDEBUG


#endif //INCLUDEGUARD_DEBUG_HPP
