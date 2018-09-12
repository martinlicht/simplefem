#ifndef INCLUDEGUARD_LOGGING_HPP
#define INCLUDEGUARD_LOGGING_HPP

#include <ostream>

// 
#define PING std::clog << "ping:" << __FILE__ << ":" << __LINE__ << std::endl;




// Returns a thing into which you can << all sorts of stuff.
// Example usage:
//     LOG << "This is a short message with a number: " << 5;      
//     ERR << "This is an error message.";      

#define LOG     clog
#define ERR     cerr




// Treat the following macros as PRINT 'str' commands
// Example usage:
//     NOTE "This is a short information"

#define NOTE    clog <<

#define WARN    cerr <<
#define ALERT   cerr <<
#define ERROR   cerr <<



inline void ping() { std::clog << "ping" << std::endl; }
inline void pong() { std::clog << "pong" << std::endl; }
inline void peng() { std::clog << "peng" << std::endl; }
inline void pang() { std::clog << "pang" << std::endl; }
inline void pung() { std::clog << "pung" << std::endl; }


static std::ostream* lognotice = &std::clog;
static std::ostream* loginfo   = &std::clog;

static std::ostream* logwarn   = &std::cerr;
static std::ostream* logalert  = &std::cerr;
static std::ostream* logerr    = &std::cerr;











#endif
