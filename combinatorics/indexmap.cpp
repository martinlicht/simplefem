
#include "indexmap.hpp"

#include <cassert>
#include <iostream>

#include "indexrange.hpp"




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
			return false;
			
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

bool IndexMap::isstrictlyascending() const
{
    for( int a = getSourceRange().getlow(); a < getSourceRange().gethigh(); a++ )
        if( values.at( a - getDestRange().getlow() ) >= values.at( a - getDestRange().getlow() + 1 ) )
            return false;
    return true;
}
                

IndexMap IndexMap::skip( int i ) const 
{
	assert( src.contains(i) );
	IndexMap ret( *this );
	ret.src = IndexRange( src.getlow(), src.gethigh() - 1 );
	ret.values.erase( ret.values.begin() + ( i - src.getlow() ) );
	return ret;
}

IndexMap IndexMap::attachbefore( int to ) const 
{
	assert( dest.contains(to) );
	IndexMap ret( *this );
	ret.src = IndexRange( src.getlow() - 1, src.gethigh() );
	ret.values.insert( ret.values.begin(), to );
}


