
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










