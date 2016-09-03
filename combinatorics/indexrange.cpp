

#include "indexrange.hpp"

#include <cassert>

#include <iostream>


IndexRange::IndexRange( int l, int h )
: minimum(l), maximum(h)
{}

void IndexRange::check() const
{
    // NOTE: Nothing to do here so far ...
}


void IndexRange::print( std::ostream& os ) const
{
    os << "Index Range: [ " << minimum << " : " << maximum << " ]" << std::endl;
}

int IndexRange::min() const
{
    return minimum;
}

int IndexRange::max() const
{
    return maximum;
}

int IndexRange::getlength() const
{
    return std::max( 0, maximum - minimum + 1 );
}

int IndexRange::cardinality() const
{
    return getlength();
}

bool IndexRange::isempty() const
{
    return minimum > maximum;
}

bool IndexRange::contains( int i ) const
{
    return minimum <= i && i <= maximum;
}

bool IndexRange::contains( const IndexRange& subir ) const
{
    return minimum <= subir.minimum && subir.maximum <= maximum;
}

bool IndexRange::operator== ( const IndexRange& other ) const
{
    return this->minimum == other.minimum && this->maximum == other.maximum;
}

int IndexRange::element2position( int i ) const
{
    assert( contains(i) );
    return i - minimum;
}

int IndexRange::position2element( int i ) const
{
    assert( 0 <= i && i <= maximum - minimum );
    int ret = i + minimum;
    assert( contains(ret) );
    return ret;
}



