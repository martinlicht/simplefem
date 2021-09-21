#ifndef INCLUDEGUARD_DEBUG_HPP
#define INCLUDEGUARD_DEBUG_HPP

void abort();

#ifdef FLAG_USE_ORIGINAL_ASSERT_MACRO

#include <cassert>
#define Assert(x) assert(x)

#else // FLAG_USE_ORIGINAL_ASSERT_MACRO

#ifdef NDEBUG
#define Assert(x) (static_cast<void>0)
#else // NDEBUG
#define Assert(x) (static_cast<bool>(x)?(void(0)):myAssert(#x,__FILE__,__LINE__))
#endif //NDEBUG

#endif //FLAG_USE_ORIGINAL_ASSERT_MACRO



#include <cstdio>

# ifdef USE_BACKTRACER
#include <stdlib.h>
#include <execinfo.h>
#endif // USE_BACKTRACER


inline void myAssert( const char* expression, const char* filename, const int linenumber )
{
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );
    fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" );
    fprintf( stderr, "!!\n" );
    fprintf( stderr, "!!\tThe following assertion failed.\n" );
    fprintf( stderr, "!!\t%s,l.%d\n", filename, linenumber );
    fprintf( stderr, "!!\n" );
    fprintf( stderr, "!!\t\t%s\n", expression );
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



#define unreachable() \
        fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), \
        fprintf( stderr, "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), \
        fprintf( stderr, "!!\n" ), \
        fprintf( stderr, "!!\tUnreachable code reached:\n!!!!\t%s:%d\n", __FILE__, __LINE__ ), \
        fprintf( stderr, "!!\n" ), \
        fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), \
        fprintf( stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ), \
        abort()



// inline void Check(){}

// template<typename T>
// inline void Check( const char* expression, const char* filename, const int linenumber )
// {
//     LOG << arg;
// }

// template<typename T, typename... Ts>
// inline void Check( T arg, Ts... args )
// {
//     LOG << arg;
//     lg( args... );
// }

// #define unreachable 
//     []() -> void{  
//         fprintf( stderr, "Unreachable code reached: %s:%d\n", __FILE__, __LINE__ );  
//         abort();  
//         }
// // __builtin_unreachable


#endif //INCLUDEGUARD_DEBUG_HPP
