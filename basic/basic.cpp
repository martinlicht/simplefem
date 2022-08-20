
#include <vector>

#include "basic.hpp"

template class std::vector<int>;
template class std::vector<Float>;


// Since all literals throughout are double unless marked otherwise 
// we enforce that `Float` is at least enough to store double.
// Any of those should do:
// 
// static_assert( Float(std::numeric_limits<double>::max()) == std::numeric_limits<double>::max(), "Float must be at least double" );
static_assert( sizeof(Float) >= sizeof(double), "Float must be at least double" );


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
