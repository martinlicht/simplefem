
#include "multiindex.hpp"

#include <vector>
#include <iostream>

#include "../basic.hpp"

#include "indexrange.hpp"



MultiIndex::MultiIndex( const IndexRange& ir )
: IndexMap( ir, NonNegativeIntegers, std::vector<int>( ir.cardinality(), 0 ) )
{
    check();
}

MultiIndex::MultiIndex( const IndexRange& ir, const std::vector<int>& vals )
: IndexMap( ir, NonNegativeIntegers, vals )
{
    assert( ir.cardinality() == vals.size() );
    check();
}


MultiIndex::MultiIndex( const IndexRange& ir, const std::function<int(int)>& generator )
: IndexMap( ir, generator )
{
    check();
}


MultiIndex::MultiIndex( const IndexRange& ir, const std::initializer_list<int>& il )
: IndexMap( ir, il )
{
    check();
}
        





void MultiIndex::check() const
{
    #ifdef NDEBUG
    return;
    #endif
    
    IndexMap::check();
}

void MultiIndex::print( std::ostream& os, bool embellish ) const
{
    check();
    for( int p : getIndexRange() )
        os << at( p ) << "\t";
    if( embellish ) 
        os << std::endl;
}

IndexRange MultiIndex::getIndexRange() const
{
    check();
    return IndexMap::getSourceRange();
}

const std::vector<int>& MultiIndex::getvalues() const
{
    check();
    return IndexMap::getvalues();
}
        





// const int& MultiIndex::at( int p ) const 
// {
//     check();
//     assert( getSourceRange().contains(p) );
//     return getvalues().at( getSourceRange().element2position(p) );
// }
// 
// int& MultiIndex::at( int p )
// {
//     check();
//     assert( getSourceRange().contains(p) );
//     return getvalues().at( getSourceRange().element2position(p) );
// }
// 
// const int& MultiIndex::operator[]( int p ) const 
// {
//     check();
//     assert( getSourceRange().contains(p) );
//     return getvalues().at( getSourceRange().element2position(p) );
// }
// 
// int& MultiIndex::operator[]( int p )
// {
//     check();
//     assert( getSourceRange().contains(p) );
//     return getvalues().at( getSourceRange().element2position(p) );
// }






int MultiIndex::absolute() const
{
    check();
    int ret = 0;
    for( int p : getIndexRange()  )
        ret += ::absolute<int>( at( p ) );
    return ret;
}

int MultiIndex::factorial() const
{
    check();
    int ret = 1;
    for( int p : getIndexRange()  )
        ret *= factorial_integer( at( p ) );
    return ret;
}






void MultiIndex::add( int p )
{
    check();
    assert( getIndexRange().contains(p) );
    at( p )++;
}

void MultiIndex::sub( int p )
{
    check();
    assert( getIndexRange().contains(p) );
    assert( at( p ) > 0 );
    at( p )--;
}



void MultiIndex::add( int p, int n )
{
    check();
    assert( getIndexRange().contains(p) );
    at( p ) += n;
}

void MultiIndex::sub( int p, int n )
{
    check();
    assert( getIndexRange().contains(p) );
    assert( at( p ) >= n );
    at( p ) -= n;
}




void MultiIndex::add( const MultiIndex& mi )
{
    check();
    mi.check();
    assert( getIndexRange() == mi.getIndexRange() );
    assert( comparablewith( mi ) );
    for( int p : getIndexRange() )
        add( p, mi.at( p ) );
    check();
}

void MultiIndex::sub( const MultiIndex& mi )
{
    check();
    mi.check();
    assert( getIndexRange() == mi.getIndexRange() );
    assert( comparablewith( mi ) );
    for( int p : getIndexRange() )
        sub( p , mi.at( p ) );
    check();
}






bool MultiIndex::comparablewith( const MultiIndex& mi ) const 
{
    check();
    mi.check();
    return getIndexRange() == mi.getIndexRange();
}


bool MultiIndex::less( const MultiIndex& mi ) const 
{
    check();
    mi.check();
    assert( comparablewith( mi ) );
    for( int p : getIndexRange() )
        if( at( p ) >= mi.at( p ) )
            return false;
    return true;
}

bool MultiIndex::equals( const MultiIndex& mi ) const
{
    check();
    mi.check();
    assert( comparablewith( mi ) );
    for( int p : getIndexRange() )
        if( at( p ) != mi.at( p ) )
            return false;
    return true;
}

