

#include <cassert>
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
    // NOTE: Nothing to do here yet.
}


void Coordinates::check() const
{
    assert( dimension >= 0 && number >= 0 );
}



void Coordinates::print( std::ostream& os ) const
{
    for( int n = 0; n < number; n++ ) {
        for( int d = 0; d < dimension; d++ )
            os << getdata( n, d );
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


Float Coordinates::getdata( int d, int n ) const
{
    assert( 0 <= n && n < number && 0 <= d && d < dimension );
    return data.at( n * dimension + d );
}

void Coordinates::setdata( int d, int n, Float v )
{
    assert( 0 <= n && n < number && 0 <= d && d < dimension );
    data.at( n * dimension + d ) = v;
}


FloatVector Coordinates::getvector( int n ) const
{
    assert( 0 <= n && n < number );
    FloatVector ret( dimension );
    for( int d = 0; d < dimension; d++ )
        ret[ d ] = data.at( n * dimension + d );
    return ret;
}

void Coordinates::setvector( int n, const FloatVector& input ) 
{
    assert( 0 <= n && n < number );
    for( int d = 0; d < dimension; d++ )
        data.at( n * dimension + d ) = input[ d ];
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
        FloatVector temp = getvector( n );
        temp += add;
        setvector( n, temp );
    }
}

void Coordinates::lineartransform( const LinearOperator& op )
{
    assert( op.getdimin() == dimension );
    for( int n = 0; n < number; n++ ) {
        FloatVector temp = getvector( n );
        temp = op * temp;
        setvector( n, temp );
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
    data.insert( data.end(), t.begin(), t.end() );
    number++;
}

