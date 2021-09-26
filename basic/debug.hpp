#ifndef INCLUDEGUARD_DEBUG_HPP
#define INCLUDEGUARD_DEBUG_HPP

// void abort();
#include <cstdlib>

#ifdef FLAG_USE_ORIGINAL_ASSERT_MACRO

#include <cassert>

#else // FLAG_USE_ORIGINAL_ASSERT_MACRO

#ifdef NDEBUG
#define Assert(x,...) (static_cast<void>0)
#else // NDEBUG
//#define Assert(x) (static_cast<bool>(x)?(void(0)):myActualAssert(#x,__FILE__,__LINE__))
//#define Assert(x,...) (static_cast<bool>(x)?(void(0)):myTemplatedAssert( #x, __FILE__, __LINE__ __VA_OPT__(,) __VA_ARGS__) )
#define Assert(x,...) (static_cast<bool>(x)?(void(0)):myActualAssert( #x, __FILE__, __LINE__, Concat2String(__VA_ARGS__) ) )
#endif //NDEBUG

#endif //FLAG_USE_ORIGINAL_ASSERT_MACRO











#include <cstdio>

# ifdef USE_BACKTRACER
#include <stdlib.h>
#include <execinfo.h>
#endif // USE_BACKTRACER


inline void myActualAssert( const char* expression, const char* filename, const int linenumber, const std::string message = "" )
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

#ifdef USE_BACKTRACER
    {
        const int limit = 10;
        void* array[limit];
        
        int size = backtrace( array, limit );
        char** strings = backtrace_symbols( array, size );

        if( strings != NULL )
        {
            printf( "Obtained %d stack frames.\n", size );
            for ( int i = 0; i < size; i++ ) printf( "%s\n", strings[i] );
        }

        free (strings);
    }
#endif // USE_BACKTRACER
    
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
std::string Concat2String()
{
    return "";
}

// one argument to stringify, base of induction 
template< typename T >
std::string Concat2String( const T& t )
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

// recursively build a string by stringifying arguments,
// with separators in between. Base case has only one argument.
template< typename T, typename... Params >
std::string Concat2String( const T& t, const Params&... params )
{
    std::stringstream ss;
    ss << t << '\t' << Concat2String( params... );
    return ss.str();
}

template< typename... Params >
inline void myTemplatedAssert( const char* expression, const char* filename, const int linenumber, const Params&... params )
{
    myActualAssert( expression, filename, linenumber, Concat2String( params... ) );
}



#define unreachable() \
        fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), \
        fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), \
        fprintf( stderr, "!!\n" ), \
        fprintf( stderr, "!!\tUnreachable code reached:\n!!!!\t%s:%d\n", __FILE__, __LINE__ ), \
        fprintf( stderr, "!!\n" ), \
        fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), \
        fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), \
        abort()

// #define unreachable 
//     []() -> void{  
//         fprintf( stderr, "Unreachable code reached: %s:%d\n", __FILE__, __LINE__ );  
//         abort();  
//         }
// // __builtin_unreachable


#endif //INCLUDEGUARD_DEBUG_HPP
