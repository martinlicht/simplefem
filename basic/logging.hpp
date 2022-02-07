#ifndef INCLUDEGUARD_LOGGING_HPP
#define INCLUDEGUARD_LOGGING_HPP

#include <iostream>
#include <string>
#include <sstream>

#include "basic.hpp"



// #include "logger.hpp"
// #include "prefixbuffer.hpp"

/* forward declarations */
std::string protocolprefixnow();


class Logger
{
    private:
        std::ostream& internalstream;
        std::string prefix;
        bool pad_newline_if_there_is_none;
        std::string filename;
        int linenumber;
    
        bool print_file_and_line = false;
        std::stringstream internalbuffer;
        
    public:
    
        explicit inline Logger( 
            std::ostream& os,
            const std::string& prefix = "",
            const bool do_newline = false,
            const char* filename = "UNKNOWN",
            const int linenumber = -1
        )
        : internalstream( os ), prefix( prefix ), pad_newline_if_there_is_none( do_newline )
        {}

        inline ~Logger()
        {
            
            const auto str = internalbuffer.str();
            
            bool use_prefix_next  = true;
            
            if( str.empty() ) {
                internalstream << prefix;
                std::cout << "\nEMPTY\n";
                return;
            }
            
            for( int c = 0; c < str.size(); c++ )
            {

                if( use_prefix_next ) { 
                    use_prefix_next = false;
                    internalstream << prefix;
                }
                
                auto character = str.at(c);
                
                internalstream << character;

                if( character == '\n' ) 
                { 
                    internalstream.flush();
                    use_prefix_next = true;
                }
                
            }
            
            if( pad_newline_if_there_is_none && ( str.empty() || str.back() != '\n' ) )
                internalstream << nl;

            if( print_file_and_line ) {
                internalstream << prefix;
                // internalstream << "\e[91m" << filename << ':' << linenumber << "\e[39m" << '\n';
                internalstream << "" << filename << ':' << linenumber << '\n';
            }

            #ifndef NDEBUG
            internalstream.flush();
            #endif
            
        }

        template<class T>
        Logger& operator<<( const T& t )
        {
            internalbuffer << t;
            return *this;
        }
        
        Logger& operator<<( std::ostream& (*const f)(std::ostream&) )
        {
            f( internalbuffer );
            return *this;
        }
        
};





template< typename L, typename... Params >
void printf_into_logger( L logger, const char* formatstring, Params... args )
{
    logger << printf_into_string( formatstring, args... );
    
    // std::size_t length = std::snprintf(nullptr, 0, formatstring, args... ) + 1;
    // char* str = new char[length];
    // std::snprintf( str, length, formatstring, args... );
    
    // logger << str;

    // delete[] str;
}







#ifndef FLAG_USE_PRIMITIVE_LOGGING
// returns a temporary logger to write stuff to, and line breaks on destruction 
// Example usage:
//     LOG << "This is a short message with a number: " << 5;      
//     ERR << "This is an error message.";      

#define LOG     Logger( std::cout, protocolprefixnow(), false, __FILE__, __LINE__ )
#define ERR     Logger( std::cerr, protocolprefixnow(), false, __FILE__, __LINE__ )

#else 

#define LOG     std::cout
#define LOG     std::cerr

#endif 



// utilize the printf template for stream-like objects 

#define LOGPRINTF( formatstring, ...) printf_into_logger( LOG, formatstring __VA_OPT__(,) __VA_ARGS__ );
#define ERRPRINTF( formatstring, ...) printf_into_logger( ERR, formatstring __VA_OPT__(,) __VA_ARGS__ );




// treat the following macros as PRINT 'str' commands
// Example usage:
//     NOTE "This is a note"
//     WARN "This is a warning"
//     ALERT "This is an alert"
//     ERROR "This is an error"

#define NOTE    Logger( std::cout, protocolprefixnow(), true, __FILE__, __LINE__ ) <<
#define NOTICE  Logger( std::cout, protocolprefixnow(), true, __FILE__, __LINE__ ) <<

#define WARNING Logger( std::cerr, protocolprefixnow(), true, __FILE__, __LINE__ ) <<
#define ALERT   Logger( std::cerr, protocolprefixnow(), true, __FILE__, __LINE__ ) <<
#define ERROR   Logger( std::cerr, protocolprefixnow(), true, __FILE__, __LINE__ ) <<



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
