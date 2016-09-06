

#include "multiindex.hpp"


#include <vector>
#include <iostream>

#include "../basic.hpp"
#include "indexrange.hpp"



MultiIndex::MultiIndex( const IndexRange& ir )
: range(ir), values( ir.getlength() )
{
    //
}
		
void MultiIndex::check() const
{
    range.check();
    for( int p = range.min(); p <= range.max(); p++ )
        attest( values[ range.element2position(p) ] >= 0 );
}

void MultiIndex::print( std::ostream& os ) const
{
    range.check();
    for( int p = range.min(); p <= range.max(); p++ )
        os << values[range.element2position(p)] << "\t";
    os << std::endl;
}

IndexRange MultiIndex::getIndexRange() const
{
    return range;
}

const int& MultiIndex::operator[]( int p ) const 
{
    attest( range.contains(p) );
    return values[range.element2position(p)];
}

int& MultiIndex::operator[]( int p )
{
    attest( range.contains(p) );
    return values[range.element2position(p)];
}

int MultiIndex::absolute() const
{
    int ret = 0;
    for( int p = range.min(); p <= range.max(); p++ )
        ret += ::absolute<int>( values[range.element2position(p)] );
    return ret;
}

int MultiIndex::factorial() const
{
    int ret = 1;
    for( int p = range.min(); p <= range.max(); p++ )
        ret *= ::factorial<int>( values[range.element2position(p)] );
    return ret;
}

void MultiIndex::operator+=( int p )
{
    attest( range.contains(p) );
    values[ range.element2position(p) ]++;
}

void MultiIndex::operator+=( const MultiIndex& mi )
{
    attest( range == mi.getIndexRange() );
    for( int p = range.min(); p <= range.max(); p++ )
        values[ range.element2position(p) ] += mi.values[ range.element2position(p) ];
}

void MultiIndex::operator-=( int p )
{
    attest( range.contains(p) );
    values[ range.element2position(p) ]--;
}

void MultiIndex::operator-=( const MultiIndex& mi )
{
    attest( range == mi.getIndexRange() );
    for( int p = range.min(); p <= range.max(); p++ )
        values[ range.element2position(p) ] -= mi.values[ range.element2position(p) ];
}
	
bool MultiIndex::smallerthan( const MultiIndex& mi ) const 
{
    for( int p = range.min(); p <= range.max(); p++ )
        if( values[ range.element2position(p) ] >= mi.values[ range.element2position(p) ] )
            return false;
    return true;
}

bool MultiIndex::equals( const MultiIndex& mi ) const
{
    for( int p = range.min(); p <= range.max(); p++ )
        if( values[ range.element2position(p) ] != mi.values[ range.element2position(p) ] )
            return false;
    return true;
}

