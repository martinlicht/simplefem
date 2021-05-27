#ifndef INCLUDEGUARD_LOGGING_HPP
#define INCLUDEGUARD_LOGGING_HPP

#include <ostream>

#include "logger.hpp"
#include "prefixbuffer.hpp"


// returns a temporary logger to write stuff to, and line breaks on destruction 
// Example usage:
//     LOG << "This is a short message with a number: " << 5;      
//     ERR << "This is an error message.";      

#define LOG     Logger( std::cerr, protocolprefixnow(), "", __FILE__, __LINE__ )
#define ERR     Logger( std::cerr, protocolprefixnow(), "", __FILE__, __LINE__ )


// treat the following macros as PRINT 'str' commands
// Example usage:
//     NOTE "This is a note"
//     WARN "This is a warning"
//     ALERT "This is an alert"
//     ERROR "This is an error"

#define NOTE    Logger( std::cout, protocolprefixnow(), "\n", __FILE__, __LINE__ ) <<
#define NOTICE  Logger( std::cout, protocolprefixnow(), "\n", __FILE__, __LINE__ ) <<

#define WARN    Logger( std::cerr, protocolprefixnow(), "\n", __FILE__, __LINE__ ) <<
#define WARNING Logger( std::cerr, protocolprefixnow(), "\n", __FILE__, __LINE__ ) <<
#define ALERT   Logger( std::cerr, protocolprefixnow(), "\n", __FILE__, __LINE__ ) <<
#define ERROR   Logger( std::cerr, protocolprefixnow(), "\n", __FILE__, __LINE__ ) <<


// emit the current file and line number into the log stream 
// Example usage:
//     PING;

#define PING LOG << "PING: " << __FILE__ << ":" << __LINE__;







////////////////////////////////////////////
// 
//      logging via variadic templates
// 
////////////////////////////////////////////

inline void lg(){}

template<typename T>
inline void lg( T arg )
{
    LOG << arg;
}

template<typename T, typename... Ts>
inline void lg( T arg, Ts... args )
{
    LOG << arg;
    lg( args... );
}




// LEGACY DEFINITIONS:

// inline void ping() { std::clog << "ping" << std::endl; }
// inline void pong() { std::clog << "pong" << std::endl; }
// inline void peng() { std::clog << "peng" << std::endl; }
// inline void pang() { std::clog << "pang" << std::endl; }
// inline void pung() { std::clog << "pung" << std::endl; }
// 
// 
// static std::ostream* lognotice = &std::clog;
// static std::ostream* loginfo   = &std::clog;
// 
// static std::ostream* logwarn   = &std::cerr;
// static std::ostream* logalert  = &std::cerr;
// static std::ostream* logerr    = &std::cerr;

#endif
