

#include "indexrange.hpp"

#include <iostream>


IndexRange::IndexRange( int l, int h )
: minimum(l), maximum(h)
{
  check();
}

void IndexRange::check() const
{
    assert( minimum > std::numeric_limits<decltype(minimum)>::min() );
    assert( minimum < std::numeric_limits<decltype(minimum)>::max() );
    assert( maximum > std::numeric_limits<decltype(maximum)>::min() );
    assert( maximum < std::numeric_limits<decltype(maximum)>::max() );
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
    attest( contains(i) );
    return i - minimum;
}

int IndexRange::position2element( int i ) const
{
    attest( 0 <= i && i <= maximum - minimum );
    int ret = i + minimum;
    attest( contains(ret) );
    return ret;
}



