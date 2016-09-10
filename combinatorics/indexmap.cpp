
#include "indexmap.hpp"

#include <iostream>

#include "indexrange.hpp"




IndexMap::IndexMap( const IndexRange& range )
: IndexMap( range, range )
{
    check();
}

IndexMap::IndexMap( const IndexRange& from, const IndexRange& to )
: src(from), dest(to), values( std::max( src.max() - src.min() + 1, 0 ), to.min() )
{
    std::cout << "Index Map initialized without actual values" << std::endl;
    src.check();
    to.check();
    check();
}

IndexMap::IndexMap( const IndexRange& range, const std::vector<int>& values )
: IndexMap( range, range, values )
{
    check();
}

IndexMap::IndexMap( const IndexRange& from, const IndexRange& to, const std::vector<int>& values )
: src(from), dest(to), values(values)
{
    check();
}

IndexMap::IndexMap( const IndexRange& range, const std::function<int(int)>& generator )
: IndexMap( range, range, generator )
{
    check();
}

IndexMap::IndexMap( const IndexRange& from, const IndexRange& to, const std::function<int(int)>& generator )
: src(from), dest(to), values( std::max( src.max() - src.min() + 1, 0 ), to.min() )
{
    for( int e = src.min(); e <= src.max(); e++ )
        values.at( src.element2position(e) ) = generator(e);
    check();
}

IndexMap::IndexMap( const IndexRange& range, const std::initializer_list<int>& values )
: IndexMap( range, std::vector<int>(values) )
{
  check();
}

IndexMap::IndexMap( const IndexRange& from, const IndexRange& to, const std::initializer_list<int>& values )
: IndexMap( from, to, std::vector<int>(values) )
{
  check();
}


        
        



void IndexMap::check() const 
{ 
    src.check();
    dest.check();
    assert( std::max( src.max() - src.min() + 1, 0 ) == values.size() );
    for( int a = src.min(); a <= src.max(); a++ )
        assert( dest.contains( values.at( a - src.min() ) ) );
}

void IndexMap::print( std::ostream& os ) const 
{
    os << "From" << std::endl << getSourceRange();
    os << "To" << std::endl << getDestRange();
    for( int i = getSourceRange().min(); i <= getSourceRange().max(); i++ )
        os << i << " -> " << values.at( i - src.min() ) << std::endl;
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




int& IndexMap::at( int i )
{
    return (*this)[i];
}

const int& IndexMap::at( int i ) const
{
    return (*this)[i];
}

int& IndexMap::operator[]( int i )
{
    check();
    assert( src.contains(i) );
    assert( 0 <= i - src.min() );
    assert( i - src.min() < values.size() );
    return values.at( i - src.min() );
}

const int& IndexMap::operator[]( int i ) const
{
    check();
    assert( src.contains(i) );
    assert( 0 <= i - src.min() );
    assert( i - src.min() < values.size() );
    return values.at( i - src.min() );
}




bool IndexMap::isempty() const 
{
    check();
    return getSourceRange().isempty();
}


bool IndexMap::isinjective() const 
{
    check();
    for( int a = src.min(); a <= src.max(); a++ )
        for( int b = src.min(); b <= src.max(); b++ )
            if( a != b && values.at( a - src.min() ) == values.at( b - src.min() ) )
                return false;

    return true;
}

bool IndexMap::issurjective() const 
{
    check();
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
    check();
    return isinjective() && issurjective();
}

bool IndexMap::isstrictlyascending() const
{
    check();
    for( int a = src.min(); a < src.max(); a++ )
        if( values.at( a - src.min() ) >= values.at( a - src.min() + 1 ) )
            return false;
    return true;
}
          

bool IndexMap::rangecontains( int p ) const
{
    check();
    assert( getDestRange().contains(p) );
    for( int i : src )
        if( at(i) == p )
            return true;
    return false;
} 

int IndexMap::rangeposition( int p ) const
{
    check();
    assert( getDestRange().contains(p) );
    for( int i : src )
        if( at(i) == p )
            return p;
    assert(false);
} 
        
bool IndexMap::comparablewith( const IndexMap& im ) const
{
    check();
    return getSourceRange() == im.getSourceRange() && getDestRange() == im.getDestRange();
} 

bool IndexMap::equals( const IndexMap& im ) const
{
    check();
    im.check();
    assert( comparablewith( im ) );
    return getvalues() == im.getvalues();
} 

bool IndexMap::less( const IndexMap& im ) const
{
    check();
    im.check();
    assert( comparablewith( im ) );
    return getvalues() < im.getvalues();
} 



// IndexMap IndexMap::skip( int i ) const 
// {
//     check();
//     assert( src.contains(i) );
//     IndexMap ret( *this );
//     ret.src = IndexRange( src.min(), src.max() - 1 );
//     ret.values.erase( ret.values.begin() + ( i - src.min() ) );
//     ret.check();
//     return ret;
// }
// 
// IndexMap IndexMap::attachbefore( int to ) const 
// {
//     check();
//     assert( dest.contains(to) );
//     IndexMap ret( *this );
//     ret.src = IndexRange( src.min() - 1, src.max() );
//     ret.values.insert( ret.values.begin(), to );
//     ret.check();
//     return ret;
// }


