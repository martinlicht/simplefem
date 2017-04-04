

#include "multiindex.hpp"


#include <vector>
#include <iostream>

#include "../basic.hpp"
#include "indexrange.hpp"



MultiIndex::MultiIndex( const IndexRange& ir )
: range(ir), values( ir.getlength() )
{
    check();
}

MultiIndex::MultiIndex( const IndexRange& ir, const std::vector<int>& vals )
: range(ir), values( vals )
{
    assert( ir.cardinality() == vals.size() );
    check();
}

void MultiIndex::check() const
{
    range.check();
    assert( range.cardinality() == values.size() );
    for( int p = range.min(); p <= range.max(); p++ )
        assert( values.at( range.element2position(p) ) >= 0 );
}

void MultiIndex::print( std::ostream& os ) const
{
    check();
    for( int p = range.min(); p <= range.max(); p++ )
        os << values.at( range.element2position(p) ) << "\t";
    os << std::endl;
}

IndexRange MultiIndex::getIndexRange() const
{
    check();
    return range;
}

const int& MultiIndex::at( int p ) const 
{
    check();
    assert( range.contains(p) );
    return values.at( range.element2position(p) );
}

int& MultiIndex::at( int p )
{
    check();
    assert( range.contains(p) );
    return values.at( range.element2position(p) );
}

const int& MultiIndex::operator[]( int p ) const 
{
    check();
    assert( range.contains(p) );
    return values.at( range.element2position(p) );
}

int& MultiIndex::operator[]( int p )
{
    check();
    assert( range.contains(p) );
    return values.at( range.element2position(p) );
}






int MultiIndex::absolute() const
{
    check();
    int ret = 0;
    for( int p = range.min(); p <= range.max(); p++ )
        ret += ::absolute<int>( values[range.element2position(p)] );
    return ret;
}

int MultiIndex::factorial() const
{
    check();
    int ret = 1;
    for( int p = range.min(); p <= range.max(); p++ )
        ret *= ::factorial<int>( values[range.element2position(p)] );
    return ret;
}






void MultiIndex::operator+=( int p )
{
    check();
    assert( range.contains(p) );
    values[ range.element2position(p) ]++;
}

void MultiIndex::operator-=( int p )
{
    check();
    assert( range.contains(p) );
    values[ range.element2position(p) ]--;
}




void MultiIndex::operator+=( const MultiIndex& mi )
{
    check();
    mi.check();
    assert( range == mi.getIndexRange() );
    assert( comparablewith( mi ) );
    for( int p = range.min(); p <= range.max(); p++ )
        values[ range.element2position(p) ] += mi.values[ range.element2position(p) ];
    check();
}

void MultiIndex::operator-=( const MultiIndex& mi )
{
    check();
    mi.check();
    assert( range == mi.getIndexRange() );
    assert( comparablewith( mi ) );
    for( int p = range.min(); p <= range.max(); p++ )
        values[ range.element2position(p) ] -= mi.values[ range.element2position(p) ];
    check();
}






bool MultiIndex::comparablewith( const MultiIndex& mi ) const 
{
    check();
    mi.check();
    return range == mi.range;
}

        
bool MultiIndex::less( const MultiIndex& mi ) const 
{
    check();
    mi.check();
    assert( comparablewith( mi ) );
    for( int p = range.min(); p <= range.max(); p++ )
        if( values[ range.element2position(p) ] >= mi.values[ range.element2position(p) ] )
            return false;
    return true;
}

bool MultiIndex::equals( const MultiIndex& mi ) const
{
    check();
    mi.check();
    assert( comparablewith( mi ) );
    for( int p = range.min(); p <= range.max(); p++ )
        if( values[ range.element2position(p) ] != mi.values[ range.element2position(p) ] )
            return false;
    return true;
}

