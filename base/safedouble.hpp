#ifndef INCLUDEGUARD_BASE_SAFEDOUBLE_HPP
#define INCLUDEGUARD_BASE_SAFEDOUBLE_HPP

#include <limits>

#include "base.hpp"

class safedouble {
    double value;

public:
    
    // Constructor from float
    explicit safedouble( float v ) : value( static_cast<double>(v) ) {}

    // Constructor from double
    explicit safedouble( double v ) : value(v) {}

    // Constructor from long double with range check
    explicit safedouble( long double v ) : value( static_cast<double>(v) )
    {
        if( std::isfinite(v) )
            assert( not( v < std::numeric_limits<double>::lowest() || v > std::numeric_limits<double>::max() ) );
    }

    // Explicit conversion operator to double
    explicit operator double() const { return value; }

};

#endif // INCLUDEGUARD_BASE_SAFEDOUBLE_HPP
