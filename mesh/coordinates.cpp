

// #include <cassert>
#include <iostream>
#include <vector>
#include <iterator>

#include "coordinates.hpp"


Coordinates::Coordinates( int dimension, int number )
: 
    dimension( dimension ), 
    number( number ), 
    data( dimension * number )
{
    check();
}

Coordinates::Coordinates( int dimension, int number, const std::vector<Float>& data )
: 
    dimension( dimension ), 
    number( number ), 
    data( data )
{
    check();
}

void Coordinates::check() const
{
    assert( dimension >= 0 && number >= 0 );
    assert( data.size() == dimension * number );
}



void Coordinates::print( std::ostream& os ) const
{
    os << "dimension: " << dimension << " - #vertices: " << number << std::endl;
	for( int n = 0; n < number; n++ ) {
        for( int d = 0; d < dimension; d++ )
            os << getdata( n, d ) << " ";
        os << std::endl;
    }
}


void Coordinates::read( std::istream& is ) 
{
    for( int n = 0; n < number; n++ ) {
        for( int d = 0; d < dimension; d++ ) {
            Float temp;
            is >> temp;
            setdata( n, d, temp );
            assert( ! is.fail() );
        }
    }
}
		




int Coordinates::getdimension() const
{
    return dimension;
}

int Coordinates::getnumber() const
{
    return number;
}

IndexRange Coordinates::getIndexRange() const
{
    return IndexRange( 0, getnumber() - 1 );
}


Float Coordinates::getdata( int n, int d ) const
{
    assert( 0 <= n && n < number && 0 <= d && d < dimension );
    return data.at( n * dimension + d );
}

void Coordinates::setdata( int n, int d, Float v )
{
    assert( 0 <= n && n < number && 0 <= d && d < dimension );
    data.at( n * dimension + d ) = v;
}


FloatVector Coordinates::getvectorclone( int n ) const
{
    assert( 0 <= n && n < number );
    return getvectorclone( n, 1. );
}

FloatVector Coordinates::getvectorclone( int n, Float scale ) const
{
    assert( 0 <= n && n < number );
    FloatVector ret( dimension );
    for( int d = 0; d < dimension; d++ )
        ret[ d ] = scale * data.at( n * dimension + d );
    return ret;
}

void Coordinates::loadvector( int n, const FloatVector& input ) 
{
    loadvector( n, input, 1. );
}

void Coordinates::loadvector( int n, const FloatVector& input, Float scale ) 
{
    assert( 0 <= n && n < number );
    assert( input.getdimension() == dimension );
    for( int d = 0; d < dimension; d++ )
        data.at( n * dimension + d ) = scale * input[ d ];
}







void Coordinates::scale( Float alpha )
{
    for( int n = 0; n < number; n++ )
        for( int d = 0; d < dimension; d++ )
            data.at( n * dimension + d ) *= alpha;
}
                                
void Coordinates::shift( const FloatVector& add )
{
    assert( add.getdimension() == dimension );
    for( int n = 0; n < number; n++ ) {
        FloatVector temp = getvectorclone( n );
        temp += add;
        loadvector( n, temp );
    }
}

void Coordinates::lineartransform( const LinearOperator& op )
{
    assert( op.getdimin() == dimension );
    assert( op.getdimout() == dimension );
    for( int n = 0; n < number; n++ ) {
        FloatVector temp = getvectorclone( n );
        temp = op * temp;
        loadvector( n, temp );
    }
}



void Coordinates::append( const Coordinates& co )
{
    assert( dimension == co.getdimension() );
    data.insert( this->data.end(), co.data.begin(), co.data.end() );
}
                
void Coordinates::append( const FloatVector& v )
{
    assert( dimension == v.getdimension() );
    const std::vector<Float>& t = v.getdata();
    data.insert( this->data.end(), t.begin(), t.end() );
    number++;
}



DenseMatrix Coordinates::getLinearPart( const IndexMap& im ) const
{
    assert( im.getDestRange() == getIndexRange() );
    
    IndexRange imsrc = im.getSourceRange();
    assert( imsrc.min() == 0 && imsrc.max() <= getdimension() );
    
    DenseMatrix ret( getdimension(), std::max(0,imsrc.max()-1) );
    assert( ret.getdimout() == getdimension() );
    assert( ret.getdimin() == std::max(0,imsrc.max()-1) );
    
    for( int p = 1; p <= imsrc.max(); p++ )
        ret.setcolumn( p-1, 
            getvectorclone( im[p] ) - getvectorclone( im[0] ) 
            );
    
    return ret;
}

FloatVector Coordinates::getShiftPart( const IndexMap& im ) const
{
    assert( im.getDestRange() == getIndexRange() );
    IndexRange imsrc = im.getSourceRange();
    assert( !(im.isempty()) && imsrc.min() == 0 && imsrc.max() <= getdimension() );
    int index = im[0];
    return getvectorclone( index );
}





