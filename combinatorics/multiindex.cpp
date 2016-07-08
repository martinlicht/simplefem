

#include "multiindex.hpp"

#include <cassert>

#include <vector>
#include <iostream>

#include "../basic.hpp"
#include "indexrange.hpp"



MultiIndex::MultiIndex( IndexRange ir )
: range(ir), values( ir.getlength() )
{
	//
}
		
void MultiIndex::check() const
{
	range.check();
	for( int p = range.min(); p <= range.max(); p++ )
		assert( values[range.place(p)] >= 0 );
}

void MultiIndex::print( std::ostream& os ) const
{
	range.check();
	for( int p = range.min(); p <= range.max(); p++ )
		os << values[range.place(p)] << "\t";
	os << std::endl;
}

IndexRange MultiIndex::getIndexRange() const
{
	return range;
}

const int& MultiIndex::operator[](int p) const 
{
	assert( range.contains(p) );
	return values[range.place(p)];
}

int& MultiIndex::operator[](int p)
{
	assert( range.contains(p) );
	return values[range.place(p)];
}

int MultiIndex::absolute() const
{
	int ret = 0;
	for( int p = range.min(); p <= range.max(); p++ )
		ret += ::absolute<int>( values[range.place(p)] );
	return ret;
}

int MultiIndex::factorial() const
{
	int ret = 1;
	for( int p = range.min(); p <= range.max(); p++ )
		ret *= ::factorial<int>( values[range.place(p)] );
	return ret;
}

void MultiIndex::operator+=( int p )
{
	assert( range.contains(p) );
	values[ range.place(p) ]++;
}

void MultiIndex::operator+=( const MultiIndex& mi )
{
	assert( range == mi.getIndexRange() );
	for( int p = range.min(); p <= range.max(); p++ )
		values[ range.place(p) ] += mi.values[ range.place(p) ];
}

void MultiIndex::operator-=( int p )
{
	assert( range.contains(p) );
	values[ range.place(p) ]--;
}

void MultiIndex::operator-=( const MultiIndex& mi )
{
	assert( range == mi.getIndexRange() );
	for( int p = range.min(); p <= range.max(); p++ )
		values[ range.place(p) ] -= mi.values[ range.place(p) ];
}
	
bool MultiIndex::smallerthan( const MultiIndex& mi ) const 
{
	for( int p = range.min(); p <= range.max(); p++ )
		if( values[ range.place(p) ] >= mi.values[ range.place(p) ] )
			return false;
	return true;
}

bool MultiIndex::equals( const MultiIndex& mi ) const
{
	for( int p = range.min(); p <= range.max(); p++ )
		if( values[ range.place(p) ] != mi.values[ range.place(p) ] )
			return false;
	return true;
}

