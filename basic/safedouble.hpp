#ifndef SAFE_DOUBLE_HPP
#define SAFE_DOUBLE_HPP

#include <limits>

#include "basic.hpp"

class safedouble {
    double value;

public:
    
    // Constructor from float
    safedouble( float v ) : value( static_cast<double>(v) ) {}

    // Constructor from double
    safedouble( double v ) : value(v) {}

    // Constructor from long double with range check
    safedouble( long double v ) : value( static_cast<double>(v) )
    {
        assert( not( v < std::numeric_limits<double>::lowest() || v > std::numeric_limits<double>::max() ) );
    }

    // Explicit conversion operator to double
    explicit operator double() const { return value; }

};

#endif // SAFE_DOUBLE_HPP
