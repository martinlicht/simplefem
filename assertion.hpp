
#ifndef __CUSTOM_ASSERTION
#define __CUSTOM_ASSERTION

#include<iostream>
#include<exception>
#include<stdexcept>

#ifdef NDEBUG
# define assert(EX)
# define warning(EX)
# define raiserror()
#else
# define warnings_are_fatal false
# define assert(EX) (void)((EX) || (__customassert (#EX, nullptr, __FILE__, __LINE__,true),0))
# define warning(EX) (void)((EX) || (__customassert (#EX, nullptr, __FILE__, __LINE__,warnings_are_fatal),0))
# define raiserror() (void)(__customassert ("Unreachable code has been reached.", nullptr, __FILE__, __LINE__,true),0))
#endif

#ifndef __cplusplus
# error Error: C++ is needed.
#endif
#if __cplusplus <= 199711L
# error This library needs at least a C++11 compliant compiler
#endif


inline void __customassert( const char* exp, const char* msg, const char* file, int line, bool shallithrow )
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
  
  if( shallithrow )
#ifdef __EXCEPTIONS
    throw std::runtime_error( "Runtime error exception thrown." );
#else
    abort();
#endif 
  else
    std::cout << "Proceed with execution." << std::endl;
  
}


#endif

