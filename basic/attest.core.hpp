
#ifndef INCLUDEGUARD_ATTEST
#define INCLUDEGUARD_ATTEST

#include<iostream>
#include<exception>
#include<stdexcept>

#ifndef __cplusplus
# error Error: C++ is needed.
#endif

#if __cplusplus <= 199711L
# error This library needs at least a C++11 compliant compiler
#endif



inline int __feecpp_attest( const char* msg, 
                            const char* exp, 
                            const char* file, 
                            int linenumber, 
                            bool shallithrow = true )
{
  
  std::cout << ( msg != nullptr ? msg : "[[no message string]]" )
            << std::endl
            << ( exp != nullptr ? exp : "[[no expression string]]" )
            << std::endl
            << "File: "
            << file
            << std::endl
            << "line: "
            << linenumber
            << std::endl;
            
  std::cout.flush();
  
  if( shallithrow ) {
    abort();
  } else {
    std::cout << "Proceed with execution." << std::endl;
  }
  
  return 0;
}

// Define macros that call the custom attest 
// and which are ALWAYS compiled.

#define warnings_are_fatal false

#define enforce_attest(EX)      ((void) (__feecpp_attest( "Error: the following attestion has failed.", #EX, __FILE__, __LINE__ ),0))
#define enforce_warning(EX)     ((void) (__feecpp_attest( "Warning: the following attestion has failed.", #EX, __FILE__, __LINE__, warnings_are_fatal ),0))
#define enforce_shout(msg)      ((void) (__feecpp_attest( msg, nullptr, __FILE__, __LINE__, warnings_are_fatal ),0))
#define enforce_raiseerror(msg) ((void) (__feecpp_attest( "Error:", msg, __FILE__, __LINE__ ),0))
#define enforce_unreachable()   ((void) (__feecpp_attest( "Error: unreachable code has been reached!", "", __FILE__, __LINE__ ),0))



#endif

