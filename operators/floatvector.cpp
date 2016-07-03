
#include <cstdlib>
#include <cmath>

#include <ostream>
#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"

FloatVector::FloatVector( int dim )
: dimension(dim), data(dim)
{
	zero();
}

FloatVector::FloatVector( const FloatVector& src )
: dimension( src.getdimension() ), data( src.getdimension() )
{
	copydatafrom( src );
}

FloatVector::FloatVector( const FloatVector& src, Float alpha )
: dimension( src.getdimension() ), data( src.getdimension() )
{
	copydatafrom( alpha, src );
}

void FloatVector::check() const 
{
	assert( dimension >= 0 );
	if( dimension == 0 ) std::cout << "WARNING: VECTOR OF DIMENSION ZERO." << std::endl;
}

void FloatVector::print( std::ostream& output ) const 
{
	for( int p = 0; p < dimension; p++ )
		output << p << ": " << getentry(p) << std::endl;
}



void FloatVector::zero() 
{
	for( int p = 0; p < dimension; p++ )
		setentry( p, 0. ); 
}

void FloatVector::random() 
{
	for( int p = 0; p < dimension; p++ )
		setentry( p, sqrt( rand() ) ); 
}

void FloatVector::scale( Float alpha ) 
{
	for( int p = 0; p < dimension; p++ )
		setentry( p, alpha * getentry( p ) ); 
}

Float FloatVector::setentry( int p, Float value )
{
	assert( 0 <= p && p < dimension );
	data.at(p) = value;
	return data.at(p);
}

Float FloatVector::getentry( int p ) const 
{
	assert( 0 <= p && p < dimension );
	return data.at(p);
}

Float& FloatVector::operator[]( int p )
{
	assert( 0 <= p && p < dimension );
	return data.at(p);
}

const Float& FloatVector::operator[]( int p ) const
{
	assert( 0 <= p && p < dimension );
	return data.at(p);
}
		

	
int FloatVector::getdimension() const 
{
	return dimension;
}



void FloatVector::copydatafrom( const FloatVector& source )
{
	copydatafrom( 1., source );
}

void FloatVector::copydatafrom( Float scaling, const FloatVector& source )
{
    assert( dimension == source.getdimension() );
	for( int p = 0; p < dimension; p++ )
		setentry( p, scaling * source.getentry( p ) ); 	
}
	
	
void FloatVector::adddatafrom( const FloatVector& source )
{
	adddatafrom( 1., source );
}

void FloatVector::adddatafrom( Float scalingsource, const FloatVector& source )
{
	adddatafrom( 1., scalingsource, source );
}
	
void FloatVector::adddatafrom( Float scalingself, Float scalingsource, const FloatVector& source )
{
    assert( dimension == source.getdimension() );
	for( int p = 0; p < dimension; p++ )
		setentry( p, scalingself * getentry( p ) + scalingsource * source.getentry( p ) ); 	
}
	

const std::vector<Float>& FloatVector::getdata() const
{
    return data;
}
	
	