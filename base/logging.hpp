#ifndef INCLUDEGUARD_LOGGING_HPP
#define INCLUDEGUARD_LOGGING_HPP

#include "base.hpp"
#include "safedouble.hpp"


// ============================================================================
// Provides a simple logging stream so that we can write 
// 
//     LOG << "Value = " << value << nl;
// 
// If the variable `USE_PRIMITIVE_LOGGING` is defined, 
// then we default to the C++ streams 
// Otherwise, we use the following custom Logger class
// 
// The custom logger relies on a global state, just like the std C++ streams 
// ============================================================================




// ============================================================================
// Unless we disable custom logging, define the logging class
// ============================================================================

#ifndef USE_PRIMITIVE_LOGGING

#include <cstdio>

#include <string>

// global boolean variable to signal whether the last write finished with a fresh line 
// extern bool log_has_a_fresh_line;

class Logger final
{
    private:
        
        std::string internal = "";
        
        bool use_cerr; //std::ostream& internalstream;
        bool pad_newline_if_there_is_none;
        std::string filename;
        int linenumber;
    
        // class-scope boolean variable to signal whether the last write finished with a fresh line 
        static bool log_has_a_fresh_line;

        // static const bool print_file_and_line = false;
        
    public:
    
        // inline
        explicit Logger( 
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
        Logger(const Logger&)            = delete;
        Logger& operator=(const Logger&) = delete;
        Logger(Logger&&)                 noexcept = delete;
        Logger& operator=(Logger&&)      noexcept = delete;
        ~Logger() noexcept;


        Logger& operator<<( char input ) {
            internal += input;
            return *this;
        }

        Logger& operator<<( const char* input ) {
            internal += input;
            return *this;
        }

        Logger& operator<<( const std::string& input ) {
            internal += input;
            return *this;
        }

        Logger& operator<<( const std::string&& input ) {
            internal += input;
            return *this;
        }

        Logger& operator<<( const void* input ) {
            char buffer[ sizeof(decltype(input)) * 2 + 10 + 1 ]; // how pointers are printed is implementation-defined 
            (void_discard)std::snprintf( buffer, sizeof(buffer), "%p", input );
            internal += buffer;
            return *this;
        }

        template <typename T>
        typename std::enable_if< std::is_integral<T>::value && std::is_signed<T>::value, Logger&>::type
        operator<<(T input) {
            char buffer[ std::numeric_limits<T>::digits10+1 + 1 + 1];
            (void_discard)std::snprintf(buffer, sizeof(buffer), "%jd", static_cast<intmax_t>(input));
            internal += buffer;
            return *this;
        }
    
        template <typename T>
        typename std::enable_if< std::is_integral<T>::value && std::is_unsigned<T>::value, Logger&>::type
        operator<<(T input) {
            char buffer[ std::numeric_limits<T>::digits10+1 + 1 + 1];
            (void_discard)std::snprintf(buffer, sizeof(buffer), "%ju", static_cast<uintmax_t>(input));
            internal += buffer;
            return *this;
        }
    
        // template <typename T, typename = decltype(std::declval<T>().text())>
        // Logger& operator<<(const T& input) {
        //     std::string text = input.text();
        //     internal += text;
        //     return *this;
        // }
    
        // Logger& operator<<( intmax_t input ) {
        //     char buffer[ sizeof(decltype(input)) * 3 + 1 + 1 ];
        //     (void_discard)std::snprintf( buffer, sizeof(buffer), "%jd", input );
        //     internal += buffer;
        //     return *this;
        // }

        // Logger& operator<<( uintmax_t input ) {
        //     char buffer[ sizeof(decltype(input)) * 3 + 1 + 1 ];
        //     (void_discard)std::snprintf( buffer, sizeof(buffer), "%ju", input );
        //     internal += buffer;
        //     return *this;
        // }

        // template <typename T>
        // typename std::enable_if<std::is_floating_point<T>::value, Logger&>::type
        // operator<<( T input ) {
        //     char buffer[ std::numeric_limits<float>::max_digits10 + std::numeric_limits<float>::max_exponent10 + 10 + 1];
        //     (void_discard)std::snprintf( buffer, sizeof(buffer), "%.10lg", (double)(safedouble)input );
        //     internal += buffer;
        //     return *this;
        // }

        Logger& operator<<( float input ) {
            return *this << (double)input;
        }

        Logger& operator<<( double input ) {
            char buffer[ std::numeric_limits<double>::max_digits10 + std::numeric_limits<double>::max_exponent10 + 10 + 1];
            (void_discard)std::snprintf( buffer, sizeof(buffer), "%.10g", input );
            internal += buffer;
            return *this;
        }

        Logger& operator<<( long double input ) {
            char buffer[ std::numeric_limits<long double>::max_digits10 + std::numeric_limits<long double>::max_exponent10 + 10 + 1];
            (void_discard)std::snprintf( buffer, sizeof(buffer), "%.10lg", (double)(safedouble)input );
            internal += buffer;
            return *this;
        }

};

#endif // USE_PRIMITIVE_LOGGING






// ============================================================================
// Define macros LOG and ERR for streaming
// ============================================================================

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

#endif // USE_PRIMITIVE_LOGGING



// ============================================================================
// Define printf-like macros for logging
// ============================================================================

#if __cplusplus < 202002L
#define LOGPRINTF(...) printf_into_stream( LOG, __VA_ARGS__ );
#define ERRPRINTF(...) printf_into_stream( ERR, __VA_ARGS__ );
#else
#define LOGPRINTF( formatstring, ...) printf_into_stream( LOG, formatstring __VA_OPT__(,) __VA_ARGS__ );
#define ERRPRINTF( formatstring, ...) printf_into_stream( ERR, formatstring __VA_OPT__(,) __VA_ARGS__ );
#endif



// ============================================================================
// The following macros work as simple PRINT "str" commands
// 
//     NOTE "This is a note"
//     WARN "This is a warning"
//     ALERT "This is an alert"
//     ERROR "This is an error"
// ============================================================================

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

#endif // USE_PRIMITIVE_LOGGING




// ============================================================================
// Emit the current file and line number into the log stream.
// Use this macro like a usual statement: PING;
// ============================================================================

#define PING LOG << "PING: " << __FILE__ << ":" << __LINE__ << nl;






#endif // INCLUDEGUARD_LOGGING_HPP
