
#include <chrono>
#include <vector>

#include "basic.hpp"
#include "constants.hpp"

template class std::vector<char>;
template class std::vector<int>;
template class std::vector<Float>;


// Since all floating-point literals throughout are double unless marked otherwise 
// we enforce that `Float` is at least enough to store double.
// Any of those should do:
// 
// static_assert( Float(std::numeric_limits<double>::max()) == std::numeric_limits<double>::max(), "Float must be at least double" );
// static_assert( sizeof(Float) >= sizeof(double), "Float must be at least double" );








int count_white_space( const std::string& str ) 
{ 
    int ret = 0;  
    for( int c = 0; c < str.size(); c++ ) if( std::isspace( str[c] ) ) ret++; 
    return ret;
} 

std::string tab_each_line( std::string str ) 
{ 
    str.insert( 0, 1, '\t' );
    for( int c = str.size(); c > 0; c-- ) {
        if( str[c-1] == '\n' )
            str.insert(c, 1, '\t');
    }
    return str;
} 

// inline std::string pad_each_line( std::string str, std::string pad )
// {
//     int c = 0;
//     for( int i = 0; i < str.length(); i++ ) if( str[i] == '\n' ) c++;
//     std::string ret = pad;
//     ret.reserve( str.length() + c * pad.length() );
//     for( int i = 0; i < str.length(); i++ )
//     {
//         ret += str[i];
//         if( str[i] == '\n' ) ret += pad.length();
    
//     return ret;
// }











#include <cstdio>
#include <cstdarg>

std::string printf_into_string( const char* formatstring, ... )
{
    
    va_list args;
    
    va_start( args, formatstring );
    std::size_t length = std::vsnprintf(nullptr, 0, formatstring, args ) + 1;
    va_end( args );
    
    char* c_str = new char[length];
    
    va_start( args, formatstring );
    std::vsnprintf( c_str, length, formatstring, args );
    va_end( args );
    
    std::string ret( c_str );
    delete[] c_str;
    return ret;
}





// TODO: Move time to cpp
static_assert( std::is_integral< decltype( std::chrono::time_point_cast< std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count() ) >::value , "Time measurement must be integral" );

timestamp gettimestamp()
{
    
    static timestamp start_timestamp = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    timestamp                    now = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    Assert( now >= start_timestamp );
    
    return now - start_timestamp;
}



// TODO: move to utility 

std::string timestamp2measurement( const timestamp& t )
{
    return to_text( static_cast<uintmax_t>(t) ) + "ms";
}

// std::string measurementnow( const timestamp& t ) // TODO Remove this line 
std::string measurementnow()
{
    return timestamp2measurement( gettimestamp() );
}



std::string timestamp2digitalcode( const timestamp& t )
{
    const int numdigits = 10;
    const int fulllength = numdigits+1;
    char digits[fulllength];
    snprintf( digits, fulllength, "%*ju", numdigits, (uintmax_t)t );
    for( int i = 0; i < numdigits; i++ ) if( digits[i] == ' ' ) digits[i] = '_';
    return std::string(digits);
}

std::string digitalcodenow()
{
    return timestamp2digitalcode( gettimestamp() );
}



std::string protocolprefixnow()
{
    // static const std::string foo = std::string("\e[36m[");
    // static const std::string bar = std::string("]\e[39m\t");
    static const std::string foo = std::string("[");
    static const std::string bar = std::string("]\t");
    return foo + digitalcodenow() + bar;
}










