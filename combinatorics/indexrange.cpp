
#include "indexrange.hpp"

#include <limits>
// #include <ostream>
#include <sstream>
#include <string>


IndexRange::IndexRange( int from, int to )
: minimum(from), maximum(to)
{
  
  if( minimum > maximum ) 
      maximum = minimum - 1;
  
  IndexRange::check();
}


void IndexRange::check() const
{
    assert( std::numeric_limits<decltype(minimum)>::min() <= minimum );
    assert( minimum <= std::numeric_limits<decltype(minimum)>::max() );
    
    assert( std::numeric_limits<decltype(maximum)>::min() <= maximum );
    assert( maximum <= std::numeric_limits<decltype(maximum)>::max() );
}

std::string IndexRange::text( bool embellish ) const
{
    std::ostringstream ss;
    
    if( embellish )
        ss << '[' << minimum << " .. " << maximum << ']';
    else
        ss << minimum << ' ' << maximum;
    
    return ss.str();
}

// void IndexRange::print( std::ostream& os, bool embellish ) const
// {
//     check();
//     os << text( embellish ) << nl;
// }


int IndexRange::min() const
{
    check();
    assert( !is_empty() );
    return minimum;
}

int IndexRange::max() const
{
    check();
    assert( !is_empty() );
    return maximum;
}

int IndexRange::cardinality() const
{
    check();
    assert( ::maximum( 0, maximum - minimum + 1 ) == maximum - minimum + 1 );
    return maximum - minimum + 1;
}

bool IndexRange::is_empty() const
{
    check();
    return minimum > maximum;
}

bool IndexRange::contains( int element ) const
{
    check();
    return minimum <= element && element <= maximum;
}

bool IndexRange::contains( const IndexRange& subir ) const
{
    check();
    subir.check();
    
    if( subir.is_empty() )
        return true;
    else if( is_empty() ) // since we know that subir is not empty 
        return false;
    else
        return minimum <= subir.minimum && subir.maximum <= maximum;
}

bool IndexRange::is_equal( const IndexRange& other ) const
{
    check();
    other.check();
    
    if( is_empty() )
        return other.is_empty();
    else if( other.is_empty() )
        return false;
    else
        return this->minimum == other.minimum && this->maximum == other.maximum;
}

int IndexRange::element2position( int element ) const
{
    check();
    assert( contains(element) );
    return element - minimum;
}

int IndexRange::position2element( int position ) const
{
    check();
    assert( !is_empty() );
    assert( 0 <= position && position <= maximum - minimum );
    int ret = position + minimum;
    assert( contains(ret) );
    return ret;
}



