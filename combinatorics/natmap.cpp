

#include <cassert>
#include <iostream>

#include "natmap.hpp"


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



IndexMap::IndexMap( IndexRange from, IndexRange to )
: src(from), dest(to), values(0)
{
	src.check();
	to.check();
	values.resize( src.gethigh() - src.getlow() + 1 );
}

void IndexMap::check() const 
{ 
	src.check();
	dest.check();
	assert( src.gethigh() - src.getlow() + 1 == values.size() );
	for( int a = src.getlow(); a <= src.gethigh(); a++ )
		assert( dest.contains( values.at( a - src.getlow() ) ) );
}

void IndexMap::print( std::ostream& os ) const 
{
	os << "From" << std::endl
	   << getSourceRange() << "To" << std::endl
	   << getDestRange() << std::endl;
	for( int i = getSourceRange().getlow(); i <= getSourceRange().gethigh(); i++ )
		os << i << " -> " << (*this)[i] << std::endl;
}

const IndexRange& IndexMap::getSourceRange() const 
{
	return src;
}
		
const IndexRange& IndexMap::getDestRange() const 
{
	return dest;
}

int& IndexMap::operator[]( int i )
{
	assert( src.contains(i) );
	assert( 0 <= i - src.getlow() );
	assert( i - src.getlow() < values.size() );
	return values.at( i - src.getlow() );
}

const int& IndexMap::operator[]( int i ) const
{
	assert( src.contains(i) );
	return values.at( i - src.getlow() );
}

bool IndexMap::isinjective() const 
{
	for( int a = src.getlow(); a <= src.gethigh(); a++ )
	for( int b = src.getlow(); b <= src.gethigh(); b++ )
		if( a != b && values.at( a - src.getlow() ) == values.at( b - src.getlow() ) )
		{
			std::cout 
			<< a - src.getlow() << space
			<< b  - src.getlow() << space
			<< values.at(a) << space << values.at(b) << std::endl;
			return false;
		}
			
	return true;
}

bool IndexMap::issurjective() const 
{
	for( int a = dest.getlow(); a <= dest.gethigh(); a++ ) 
	{
		bool flag = false;
		for( int b = src.getlow(); b <= src.gethigh(); b++ )
			if( values.at( b - src.getlow() ) == a )
				flag = true;
		if( !flag )
			return false;
	}
	return true;
}

bool IndexMap::isbijective() const 
{
	return isinjective() && issurjective();
}

