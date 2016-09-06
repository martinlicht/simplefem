

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
    check();
    os << "Index Range: [ " << minimum << " : " << maximum << " ]" << std::endl;
}

int IndexRange::min() const
{
    check();
    return minimum;
}

int IndexRange::max() const
{
    check();
    return maximum;
}

int IndexRange::getlength() const
{
    check();
    return std::max( 0, maximum - minimum + 1 );
}

int IndexRange::cardinality() const
{
    check();
    return getlength();
}

bool IndexRange::isempty() const
{
    check();
    return minimum > maximum;
}

bool IndexRange::contains( int i ) const
{
    check();
    return minimum <= i && i <= maximum;
}

bool IndexRange::contains( const IndexRange& subir ) const
{
    check();
    subir.check();
    return minimum <= subir.minimum && subir.maximum <= maximum;
}

bool IndexRange::operator== ( const IndexRange& other ) const
{
    check();
    other.check();
    return this->minimum == other.minimum && this->maximum == other.maximum;
}

int IndexRange::element2position( int i ) const
{
    check();
    attest( contains(i) );
    return i - minimum;
}

int IndexRange::position2element( int i ) const
{
    check();
    attest( 0 <= i && i <= maximum - minimum );
    int ret = i + minimum;
    attest( contains(ret) );
    return ret;
}



