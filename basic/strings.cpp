
#include "strings.hpp"


#ifdef FLAG_USE_CUSTOM_STRINGS
#else // FLAG_USE_CUSTOM_STRINGS
#endif // FLAG_USE_CUSTOM_STRINGS



template std::string to_text( signed char );
template std::string to_text( signed short );
template std::string to_text( signed int );
template std::string to_text( signed long );
template std::string to_text( signed long long );
    
template std::string to_text( unsigned char );
template std::string to_text( unsigned short );
template std::string to_text( unsigned int );
template std::string to_text( unsigned long );
template std::string to_text( unsigned long long );

template std::string to_text( float );
template std::string to_text( double );
template std::string to_text( long double );







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

