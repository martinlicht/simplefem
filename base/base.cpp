
#include <cstdio>
#include <cstdarg>
#include <chrono>
#include <cctype>

#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "base.hpp"
#include "constants.hpp"

template class std::vector<char>;
template class std::vector<int>;
template class std::vector<std::size_t>;
template class std::vector<Float>;


// Since all floating-point literals throughout are double unless marked otherwise 
// we enforce that `Float` is at least enough to store double. Any of those should do:
// 
// static_assert( Float(std::numeric_limits<double>::max()) == std::numeric_limits<double>::max(), "Float must be at least double" );
// static_assert( sizeof(Float) >= sizeof(double), "Float must be at least double" );








int count_white_space( const std::string& str ) 
{ 
    int ret = 0;  
    for( int c = 0; c < str.size(); c++ ) if( 0 != std::isspace( str[c] ) ) ret++; 
    return ret;
} 

std::string tab_each_line( std::string str ) 
{ 
    str.insert( 0, 1, '\t' );
    for( int c = str.size(); c > 0; c-- ) {
        if( str[c-1] == '\n' )
            str.insert(c, 1, '\t' );
    }
    return str;
} 


int string_to_integer( const char* s, const char* __restrict__ *endptr, unsigned int base, bool& has_overflown ) {
    
    assert( s != nullptr );
    assert( ( 2 <= base and base <= 36 ) or ( base == 0 ) );

    unsigned int result = 0;
    
    int sign = 1;

    has_overflown = false;
    
    // Skip leading whitespace
    while( std::isspace( (unsigned char)*s ) ) { s++; }

    if( *s == '\0' ) {
        // printf("Only white space\n");
        return 0;
    } 
        
    // Handle optional sign
    if( *s == '-' ) { 
        sign = -1;
        s++;
    } else if( *s == '+' ) {
        s++;
    }

    if( *s == '\0' ) { 
        // printf("Only signs\n");
        return 0;
    }
        
    // Auto-detect if base equals zero 
    if( base == 0 ) {
        
        // printf("Auto-detect base\n");
        
        base = 10;

        if( *s == '0' ) { 
            
            base = 8;
            s++;
        
            if( *s == '\0' ) return 0;
            
            if( *s == 'x' || *s == 'X' ) {
                base = 16; 
                s++;
            }

        }
    }

    // Convert digits
    while( *s ) {
        
        int digit;
        
        // Find digits within the correct range 
        if( *s >= '0' && *s <= '9' ) {
            digit = *s - '0';
        } else if( *s >= 'a' && *s <= 'z' ) {
            digit = *s - 'a' + 10;
        } else if( *s >= 'A' && *s <= 'Z' ) {
            digit = *s - 'A' + 10;
        } else {
            break;
        }

        // printf("Digit detected: %u\n", digit );

        if( digit >= base ) break;

        if( not has_overflown and result > ( std::numeric_limits<int>::max() - digit) / base) {
            // printf("Overflowing... %u %u\n", result, digit );
            has_overflown = true;
        } else {
            result = result * base + (unsigned int)digit;
        }

        s++;
    }

    // Store pointer one past the last digit 
    if( endptr != nullptr ) 
    {
        *endptr = s;
    }

    // Assuming two's complement (C++20 and C23)
    if( sign == -1 and result > 1u + (unsigned int)(std::numeric_limits<int>::max()) ) has_overflown = true;
    
    if( has_overflown ) 
    {
        // printf("Overflown\n");
        return ( sign == 1 ) ? ( std::numeric_limits<int>::max() ) : ( std::numeric_limits<int>::min() );
    }

    return sign * result;
}














































std::string printf_into_string( const char* formatstring, ... )
{
    
    va_list args;
    
    va_start( args, formatstring );
    const std::size_t length = std::vsnprintf(nullptr, 0, formatstring, args ) + 1;
    va_end( args );
    
    char* c_str = new char[length];
    
    va_start( args, formatstring );
    (void_discard)std::vsnprintf( c_str, length, formatstring, args );
    va_end( args );
    
    std::string ret( c_str );
    delete[] c_str;
    return ret;
}

// template< typename... Params >
// inline std::string printf_into_string( const char* formatstring, Params... args )
// {
//     std::size_t length = std::snprintf(nullptr, 0, formatstring, args... ) + 1;
//     char* c_str = new char[length];
//     std::snprintf( c_str, length, formatstring, args... );
//     std::string ret( c_str );
//     delete[] c_str;
//     return ret;
// }





// TODO(martinlicht): simplify the time stamp interface and move it to logging, even with code duplication.

static_assert( std::is_integral< decltype( std::chrono::time_point_cast< std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count() ) >::value , "Time measurement must be integral" );

timestamp timestampnow()
{
    
    static const timestamp start_timestamp = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    const timestamp                    now = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    Assert( now >= start_timestamp );
    
    return now - start_timestamp;
}



std::string timestamp2measurement( const timestamp& t )
{
    return std::to_string( static_cast<uintmax_t>(t) ) + "ms";
}

std::string measurementnow()
{
    return timestamp2measurement( timestampnow() );
}



std::string timestamp2digitalcode( const timestamp& t )
{
    const int numdigits = 10;
    const int fulllength = numdigits+1;
    char digits[fulllength];
    assert( t < 9999999999 ); // ca. 115 days 
    (void_discard)snprintf( digits, fulllength, "%*ju", numdigits, (uintmax_t)t );
    for( int i = 0; i < numdigits; i++ ) if( digits[i] == ' ' ) digits[i] = '_';
    return std::string(digits);
}

std::string digitalcodenow()
{
    return timestamp2digitalcode( timestampnow() );
}



// std::string protocolprefixnow()
// {
//     // static const std::string foo = std::string("\e[36m[");
//     // static const std::string bar = std::string("]\e[39m\t");
//     static const std::string foo = std::string("[");
//     static const std::string bar = std::string("]\t");
//     return foo + digitalcodenow() + bar;
// }










