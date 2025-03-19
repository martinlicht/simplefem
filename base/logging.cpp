#include "logging.hpp"

#ifndef USE_PRIMITIVE_LOGGING

#include <cmath>
#include <cstdio>
#include <ctime>


#include <chrono>
#include <string>




// ============================================================================
// Functionality to obtain the current time as a string,
// measured in milliseconds and with leading space filled by being '_'
// ============================================================================

static_assert( std::is_integral< decltype( std::chrono::time_point_cast< std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count() ) >::value , 
                "Time measurement must be integral" );

static std::string get_current_time_as_string()
{
    static const auto start = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    const auto          now = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    auto t = now - start;

    Assert( t >= 0 );
    Assert( t < 99999'99999 ); // ca. 115 days 

    int numdigits = 10;
    
    char digits[numdigits+1];
    digits[numdigits] = '\0';

    numdigits--;
    do { digits[numdigits--] = '0' + ( t % 10 ); t /= 10; } while( t > 0 );
    while( numdigits >= 0 ) digits[numdigits--] = '_';
    
    return std::string(digits);
}





// ============================================================================
// Char codes for colored console output 
// ============================================================================

enum class TextColors : unsigned char
{
    black          = 30,
    black_light    = 90,
    
    red            = 31,
    red_light      = 91,
    
    green          = 32,
    green_light    = 92, 
    
    yellow         = 33, 
    yellow_bright  = 93, 
    
    blue           = 34,
    blue_bright    = 94,
    
    magenta        = 35,
    magenta_bright = 95,
    
    cyan           = 36,
    cyan_light     = 96, 
    
    white          = 37,
    white_light    = 97 
};




// ============================================================================
// Definitions for the custom logging class
// ============================================================================

bool Logger::log_has_a_fresh_line = true;

Logger::Logger( 
            bool use_cerr, //std::ostream& os,
            const bool do_newline,
            const char* filename,
            const int linenumber
        )
        : 
        use_cerr( use_cerr ), // internalstream( use_cerr ? std::cerr : std::cout ), // internalstream( os ),
        pad_newline_if_there_is_none( do_newline ),
        filename( filename ),
        linenumber( linenumber )
        {
            // *this << std::setprecision(10);
        }


// auto digitalcodenow() -> std::string;

Logger::~Logger() noexcept
{
    
    FILE* f = ( use_cerr ? stderr : stdout );
    
    const auto& message_string = internal;
    
    /* 1. complete the prefix string */

    const std::string time_code = get_current_time_as_string(); // digitalcodenow(); 

    #ifdef USE_COLORED_OUTPUT
    const std::string colorcode_begin = "\033[96m";
    const std::string colorcode_close = "\033[m";
    #else
    const std::string colorcode_begin = "";
    const std::string colorcode_close = "";
    #endif
    const std::string prefix = printf_into_string(
        "%s[%s %s %4d]%s\t", 
        colorcode_begin.c_str(), time_code.c_str(), filename.c_str(), linenumber, colorcode_close.c_str() );
    
    /* 2. introduce final output string, reserve enough spaces */

    int number_of_newlines = 0;
    for( const char& character : message_string ) 
        if( character == '\n' )
            number_of_newlines++;

    std::string output_string;

    output_string.reserve( 1 + message_string.size() + number_of_newlines * prefix.size() );

    /* 3. if we have inherited a fresh line, then insert a prefix */
    
    if( log_has_a_fresh_line ){
        log_has_a_fresh_line = false;
        output_string.insert(0,prefix);
    }
    
    /* 4. for each nl (except possibly at the final position) insert a prefix after it */
    
    for( int c = 0; c < message_string.length(); c++ ) {
    
        const char character = message_string[c];

        output_string += character;
        // fputc( character, f );

        if( character == '\n' and c != message_string.length()-1 ) 
            output_string += prefix;
        
    }
    
    /* 5. add newline if the last character is not a newline and we automatically append newlines */
    
    if( message_string.length() > 0 and message_string.back() != nl and pad_newline_if_there_is_none )
        output_string += nl;

    /* 6. memorize whether the last character is a newline */

    log_has_a_fresh_line = ( output_string.back() == nl );
    
    /* 7. output the string */

    fputs( output_string.c_str(), f );
    fflush(f);

}




#endif // USE_PRIMITIVE_LOGGING
