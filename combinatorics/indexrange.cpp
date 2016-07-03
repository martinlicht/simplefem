

#include "indexrange.hpp"

#include <cassert>

#include <iostream>


IndexRange::IndexRange( int l, int h )
: low(l), high(h)
{}

void IndexRange::check() const 
{ assert( low <= high); }

void IndexRange::print( std::ostream& os ) const 
{
	os << "Index Range: [ " << low << " : " << high << " ]" << std::endl;
}

int IndexRange::getlow() const
{
	return low;
}

int IndexRange::gethigh() const
{
	return high;
}

int IndexRange::getlength() const
{
	return high - low + 1;
}
		

bool IndexRange::contains( int i ) const 
{
	return low <= i && i <= high;
}

bool IndexRange::contains( const IndexRange& subir ) const 
{
	return low <= subir.low && subir.high <= high;
}

bool IndexRange::operator== ( const IndexRange& other ) const
{
	return this->low == other.low && this->high == other.high;
}

int IndexRange::place( int i ) const 
{
	assert( contains(i) );
	return i - low;
}




