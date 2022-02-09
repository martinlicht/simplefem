
#include <vector>

#include "basic.hpp"

template class std::vector<int>;
template class std::vector<Float>;

#include <cstdarg>

std::string printf_into_string( const char* formatstring, ... )
{
    
    va_list args;
    va_start( args, formatstring );
    std::size_t length = std::vsnprintf(nullptr, 0, formatstring, args ) + 1;
    char* c_str = new char[length];
    std::vsnprintf( c_str, length, formatstring, args );
    std::string ret( c_str );
    delete[] c_str;
    va_end( args );
    return ret;
}
