
#include <cstdlib>
#include <cmath>

#include <ostream>
#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"

FloatVector::FloatVector( int dim )
:
 data(dim)
{
	check();
}

FloatVector::FloatVector( int dim, Float initialvalue )
:
 data(dim,initialvalue)
{
    check();
}

FloatVector::FloatVector( const FloatVector& src )
: 
 data( src.getdimension() )
{
	copydatafrom( src );
    check();
}

FloatVector::FloatVector( const FloatVector& src, Float alpha )
: 
 data( src.getdimension() )
{
	copydatafrom( alpha, src );
    check();
}

FloatVector::FloatVector( int dimension, const std::function<Float(int)>& generator )
: 
 data( dimension )
{
	generatedatafrom( generator );
    check();
}

FloatVector::FloatVector( int dimension, const std::function<Float(int)>& generator, Float alpha )
: 
 data( dimension )
{
	generatedatafrom( alpha, generator );
    check();
}



void FloatVector::check() const 
{
	assert( getdimension() >= 0 );
	if( getdimension() == 0 )
        std::cout << "WARNING: VECTOR OF DIMENSION ZERO." << std::endl;
	assert( getdimension() == data.size() );
}

void FloatVector::print( std::ostream& output ) const 
{
	output << "float vector of dimension: " << getdimension() << std::endl;
    for( int p = 0; p < getdimension(); p++ )
		output << p << ": " << getentry(p) << std::endl;
}

int FloatVector::getdimension() const 
{
	return data.size();
}

Float FloatVector::setentry( int p, Float value )
{
	assert( 0 <= p && p < data.size() );
	data.at(p) = value;
	return data.at(p);
}

Float FloatVector::getentry( int p ) const 
{
	assert( 0 <= p && p < data.size() );
	return data.at(p);
}

Float& FloatVector::at( int p )
{
    return (*this)[p];
}

const Float& FloatVector::at( int p ) const
{
    return (*this)[p];
}

Float& FloatVector::operator[]( int p )
{
	assert( 0 <= p && p < data.size() );
	return data.at(p);
}

const Float& FloatVector::operator[]( int p ) const
{
	assert( 0 <= p && p < data.size() );
	return data.at(p);
}
		
const std::vector<Float>& FloatVector::getdata() const
{
    return data;
}

	




void FloatVector::zero() 
{
	for( int p = 0; p < getdimension(); p++ )
		setentry( p, 0. ); 
}

void FloatVector::random() 
{
	for( int p = 0; p < getdimension(); p++ )
		setentry( p, sqrt( rand() ) ); 
}

void FloatVector::scale( Float alpha ) 
{
	for( int p = 0; p < getdimension(); p++ )
		setentry( p, alpha * getentry( p ) ); 
}

void FloatVector::copydatafrom( const FloatVector& source )
{
	copydatafrom( 1., source );
}

void FloatVector::copydatafrom( Float scaling, const FloatVector& source )
{
    assert( getdimension() == source.getdimension() );
	for( int p = 0; p < getdimension(); p++ )
		setentry( p, scaling * source.getentry( p ) ); 	
}
	
	
void FloatVector::generatedatafrom( const std::function<Float(int)>& generator )
{
	generatedatafrom( 1., generator );
}

void FloatVector::generatedatafrom( Float scaling, const std::function<Float(int)>& generator )
{
    for( int p = 0; p < getdimension(); p++ )
		setentry( p, scaling * generator( p ) ); 	
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
    assert( getdimension() == source.getdimension() );
	for( int p = 0; p < getdimension(); p++ )
		setentry( p, scalingself * getentry( p ) + scalingsource * source.getentry( p ) ); 	
}


Float FloatVector::scalarproductwith( const FloatVector& right ) const
{
    assert( getdimension() == right.getdimension() );
    Float ret = 0.;
    for( int p = 0; p < getdimension(); p++ )
        ret += getentry(p) * right.getentry(p);
    return ret;
}

Float FloatVector::norm() const 
{
    return sqrt( scalarproductwith( *this ) );
}

Float FloatVector::maxnorm() const
{
    assert( getdimension() > 0 );
    Float ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret = std::max( ret, absolute( data.at(d) ) );
    return ret;
}

Float FloatVector::lpnorm( Float p ) const
{
    assert( p > 0 );
    assert( getdimension() > 0 );
    
    Float ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret += pow( absolute( data.at(d) ), p );
    return pow( ret, 1./p );
}





