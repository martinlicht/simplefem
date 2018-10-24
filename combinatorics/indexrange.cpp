

#include "indexrange.hpp"

#include <iostream>


IndexRange::IndexRange( int from, int to )
: minimum(from), maximum(to)
{
  if( minimum > maximum ) maximum = minimum - 1;
  check();
}

void IndexRange::check() const
{
    assert( std::numeric_limits<decltype(minimum)>::min() < minimum );
    assert( minimum < std::numeric_limits<decltype(minimum)>::max() );
    
    assert( std::numeric_limits<decltype(maximum)>::min() < maximum );
    assert( maximum < std::numeric_limits<decltype(maximum)>::max() );
}


void IndexRange::print( std::ostream& os, bool embellish ) const
{
    check();
    if( embellish )
        os << "[" << minimum << ":" << maximum << "]" << std::endl;
    else
        os << minimum << " " << maximum;
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

int IndexRange::cardinality() const
{
    check();
    return std::max( 0, maximum - minimum + 1 );
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

bool IndexRange::compare( const IndexRange& other ) const
{
    check();
    other.check();
    return this->minimum == other.minimum && this->maximum == other.maximum;
}

int IndexRange::element2position( int i ) const
{
    check();
    assert( contains(i) );
    return i - minimum;
}

int IndexRange::position2element( int i ) const
{
    check();
    assert( 0 <= i && i <= maximum - minimum );
    int ret = i + minimum;
    assert( contains(ret) );
    return ret;
}



