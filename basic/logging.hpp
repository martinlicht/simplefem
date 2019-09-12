#ifndef INCLUDEGUARD_LOGGING_HPP
#define INCLUDEGUARD_LOGGING_HPP

#include <ostream>

#include "logger.hpp"
#include "prefixbuffer.hpp"






// returns a temporary logger to write stuff to, and line breaks on destruction 
// Example usage:
//     LOG << "This is a short message with a number: " << 5;      
//     ERR << "This is an error message.";      

#define LOG     Logger( std::cout, "", "\n" )
#define ERR     Logger( std::cerr, "", "\n" )




// treat the following macros as PRINT 'str' commands
// Example usage:
//     NOTE "This is a short information"
//     WARN "This is a short information"
//     ALERT "This is a short information"
//     ERROR "This is a short information"

#define NOTE    Logger( std::cout, "", "\n" ) <<
#define NOTICE  Logger( std::cout, "", "\n" ) <<

#define WARN    Logger( std::cerr, "", "\n" ) <<
#define WARNING Logger( std::cerr, "", "\n" ) <<
#define ALERT   Logger( std::cerr, "", "\n" ) <<
#define ERROR   Logger( std::cerr, "", "\n" ) <<




// emit the current file and line number into the log stream 
// Example usage:
//     PING;

#define PING LOG << "PING: " << __FILE__ << ":" << __LINE__;






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
