
#include "indexmap.hpp"

#include <cassert>
#include <iostream>

#include "indexrange.hpp"




IndexMap::IndexMap( IndexRange from, IndexRange to )
: src(from), dest(to), values(0)
{
	src.check();
	to.check();
	values.resize( src.max() - src.min() + 1 );
}

void IndexMap::check() const 
{ 
	src.check();
	dest.check();
	assert( src.max() - src.min() + 1 == values.size() );
	for( int a = src.min(); a <= src.max(); a++ )
		assert( dest.contains( values.at( a - src.min() ) ) );
}

void IndexMap::print( std::ostream& os ) const 
{
	os << "From" << std::endl
	   << getSourceRange() << "To" << std::endl
	   << getDestRange() << std::endl;
	for( int i = getSourceRange().min(); i <= getSourceRange().max(); i++ )
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

const std::vector<int>& IndexMap::getvalues() const 
{
	return values;
}

int& IndexMap::operator[]( int i )
{
	assert( src.contains(i) );
	assert( 0 <= i - src.min() );
	assert( i - src.min() < values.size() );
	return values.at( i - src.min() );
}

const int& IndexMap::operator[]( int i ) const
{
	assert( src.contains(i) );
	return values.at( i - src.min() );
}

bool IndexMap::isinjective() const 
{
	for( int a = src.min(); a <= src.max(); a++ )
	for( int b = src.min(); b <= src.max(); b++ )
		if( a != b && values.at( a - src.min() ) == values.at( b - src.min() ) )
			return false;
			
	return true;
}

bool IndexMap::issurjective() const 
{
	for( int a = dest.min(); a <= dest.max(); a++ ) 
	{
		bool flag = false;
		for( int b = src.min(); b <= src.max(); b++ )
			if( values.at( b - src.min() ) == a )
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
    for( int a = getSourceRange().min(); a < getSourceRange().max(); a++ )
        if( values.at( a - getDestRange().min() ) >= values.at( a - getDestRange().min() + 1 ) )
            return false;
    return true;
}
          

bool IndexMap::comparablewith( const IndexMap& im ) const
{
    return getSourceRange() == im.getSourceRange() && getDestRange() == im.getDestRange();
} 

bool IndexMap::equals( const IndexMap& im ) const
{
    assert( comparablewith( im ) );
    return getvalues() == im.getvalues();
} 

bool IndexMap::less( const IndexMap& im ) const
{
    assert( comparablewith( im ) );
    return getvalues() < im.getvalues();
} 



IndexMap IndexMap::skip( int i ) const 
{
	assert( src.contains(i) );
	IndexMap ret( *this );
	ret.src = IndexRange( src.min(), src.max() - 1 );
	ret.values.erase( ret.values.begin() + ( i - src.min() ) );
	return ret;
}

IndexMap IndexMap::attachbefore( int to ) const 
{
	assert( dest.contains(to) );
	IndexMap ret( *this );
	ret.src = IndexRange( src.min() - 1, src.max() );
	ret.values.insert( ret.values.begin(), to );
}


