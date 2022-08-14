#ifndef INCLUDEGUARD_LOGGING_HPP
#define INCLUDEGUARD_LOGGING_HPP

#include "basic.hpp"



#ifndef USE_PRIMITIVE_LOGGING

#include <ostream>
#include <string>
#include <sstream>

// #include "logger.hpp"
// #include "prefixbuffer.hpp"

/* forward declarations */
std::string protocolprefixnow();


// This variable has an instance in every translation unit 
// It is not global for the entire program 
extern bool log_has_a_fresh_line;

class Logger : public std::ostringstream
{
    private:
        bool use_cerr; //std::ostream& internalstream;
        bool pad_newline_if_there_is_none;
        std::string filename;
        int linenumber;
    
        bool print_file_and_line = false;
        
    public:
    
        explicit 
        // inline
        Logger( 
            bool use_cerr, //std::ostream& os,
            const bool do_newline = false,
            const char* filename = "UNKNOWN",
            const int linenumber = -1
        )
        ;
        // : 
        // internalstream( os ),
        // pad_newline_if_there_is_none( do_newline ),
        // filename( filename ),
        // linenumber( linenumber )
        // {}

        ~Logger();

};






// class Logger2
// {
//     private:
//         std::ostream& internalstream;
//         std::string prefix;
//         bool pad_newline_if_there_is_none;
//         std::string filename;
//         int linenumber;
    
//         bool print_file_and_line = false;
//         std::ostringstream internalbuffer;
        
//     public:
    
//         explicit inline Logger2( 
//             std::ostream& os,
//             const std::string& prefix = "",
//             const bool do_newline = false,
//             const char* filename = "UNKNOWN",
//             const int linenumber = -1
//         )
//         : internalstream( os ), prefix( prefix ), pad_newline_if_there_is_none( do_newline )
//         {}

//         inline ~Logger2()
//         {
            
//             const auto str = internalbuffer.str();
//             ................. internalstream << str;
            
//         }

//         template<class T>
//         Logger& operator<<( const T& t )
//         {
//             internalbuffer << t;
//             return *this;
//         }
        
//         Logger& operator<<( std::ostream& (*const f)(std::ostream&) )
//         {
//             f( internalbuffer );
//             return *this;
//         }
        
// };





// template< typename L, typename... Params >
// void printf_into_logger( L logger, const char* formatstring, Params... args )
// {
//     logger << printf_into_string( formatstring, args... );
    
//     // std::size_t length = std::snprintf(nullptr, 0, formatstring, args... ) + 1;
//     // char* str = new char[length];
//     // std::snprintf( str, length, formatstring, args... );
    
//     // logger << str;

//     // delete[] str;
// }
#endif 






#ifndef USE_PRIMITIVE_LOGGING
// returns a temporary logger to write stuff to, and line breaks on destruction 
// Example usage:
//     LOG << "This is a short message with a number: " << 5;      
//     ERR << "This is an error message.";      

#define LOG     Logger( false, false, __FILE__, __LINE__ )
#define ERR     Logger( true, false, __FILE__, __LINE__ )

#else 

#include <iostream>

#define LOG     std::cout
#define ERR     std::cerr

#endif 



// utilize the printf template for stream-like objects 

#define LOGPRINTF( formatstring, ...) printf_into_stream( LOG, formatstring __VA_OPT__(,) __VA_ARGS__ );
#define ERRPRINTF( formatstring, ...) printf_into_stream( ERR, formatstring __VA_OPT__(,) __VA_ARGS__ );




// treat the following macros as PRINT 'str' commands
// Example usage:
//     NOTE "This is a note"
//     WARN "This is a warning"
//     ALERT "This is an alert"
//     ERROR "This is an error"

#ifndef USE_PRIMITIVE_LOGGING

#define NOTE    Logger( false, true, __FILE__, __LINE__ ) <<
#define NOTICE  Logger( false, true, __FILE__, __LINE__ ) <<

#define WARNING Logger( true, true, __FILE__, __LINE__ ) <<
#define ALERT   Logger( true, true, __FILE__, __LINE__ ) <<
#define ERROR   Logger( true, true, __FILE__, __LINE__ ) <<

#else 

#define NOTE    std::cout <<
#define NOTICE  std::cout <<

#define WARNING std::cerr <<
#define ALERT   std::cerr <<
#define ERROR   std::cerr <<

#endif 


// emit the current file and line number into the log stream 
// Example usage:
//     PING;

#define PING LOG << "PING: " << __FILE__ << ":" << __LINE__ << nl;







////////////////////////////////////////////
// 
//      logging via variadic templates
// 
////////////////////////////////////////////

inline void lg(){}

// template<typename T>
// inline void lg( T arg )
// {
//     LOG << arg << nl;
// }

template<typename T, typename... Ts>
inline void lg( const T arg, const Ts... args )
{
    LOG << arg << nl;
    lg( args... );
}







////////////////////////////////////////////
// 
//      logging via variadic templates
// 
////////////////////////////////////////////

struct OpenMP_Reporter
{
    OpenMP_Reporter();
    ~OpenMP_Reporter();
};

extern OpenMP_Reporter  omp_reporter;








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
