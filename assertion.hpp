
#ifndef __CUSTOM_ASSERTION
#define __CUSTOM_ASSERTION

#include<cassert>

#ifdef attest
#error Attest macro already defined 
#else
#define attest(EX) assert(EX)
#endif

#include<iostream>
#include<exception>
#include<stdexcept>

#ifndef __cplusplus
# error Error: C++ is needed.
#endif

#if __cplusplus <= 199711L
# error This library needs at least a C++11 compliant compiler
#endif


#ifdef NDEBUG
# define attest(EX) ((void)0)
# define warning(EX) ((void)0)
# define raiserror() ((void)0)
#else
# define warnings_are_fatal false
// # define attest(EX) ((void)((EX) || (__customassert (#EX, nullptr, __FILE__, __LINE__,true),0)))
// # define attest(EX) ( (EX) ? 0 : (__customassert (#EX, nullptr, __FILE__, __LINE__,true)) )
# define attest(EX) assert(EX)
# define warning(EX) ((void)((EX) || (__customassert (#EX, nullptr, __FILE__, __LINE__,warnings_are_fatal),0)))
# define raiseerror(msg) ((void)( __customassert ("Unreachable code has been reached.", msg, __FILE__, __LINE__,true),0))
#endif

inline int __customassert( const char* exp, const char* msg, const char* file, int line, bool shallithrow )
{
  
  if( shallithrow )
    std::cout << "ERROR: ";
  else
    std::cout << "WARNING: ";
  std::cout << file << " -- l." << line << std::endl;
  
  if( exp == nullptr )
    std::cout << "*** Empty expression string ***";
  else 
    std::cout << exp << std::endl;
  
  if( msg != nullptr )
    std::cout << msg << std::endl;
  else
    (void)0;
    
  if( shallithrow ) {
    std::cout << "Terminate execution." << std::endl;
    std::cout.flush();
// #ifdef __EXCEPTIONS
//     throw 0; 
//     std::runtime_error( "Runtime error exception thrown." );
// #else
    abort();
// #endif
  } else {
    std::cout << "Proceed with execution." << std::endl;
  }
  
  return 0;
}


#endif

