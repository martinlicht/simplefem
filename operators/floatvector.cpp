
#include <cstdlib>
#include <cmath>

#include <ostream>
#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"

FloatVector::FloatVector( int dim )
: data(dim)
{
    check();
}

FloatVector::FloatVector( int dim, Float initialvalue )
: data(dim,initialvalue)
{
    check();
}

FloatVector::FloatVector( const FloatVector& src )
: data( src.getdimension() )
{
    copydatafrom( src );
    check();
}

FloatVector::FloatVector( const FloatVector& src, Float alpha )
: data( src.getdimension() )
{
    copydatafrom( alpha, src );
    check();
}

FloatVector::FloatVector( const std::vector<Float>& vals )
: FloatVector( vals, 1. )
{
  check();
}

FloatVector::FloatVector( const std::vector<Float>& vals, Float alpha )
: FloatVector( vals.size(), [&vals](int i) -> Float{ return vals.at(i); }, alpha )
{
  check();
}

FloatVector::FloatVector( int dimension, const std::function<Float(int)>& generator )
: FloatVector( dimension, generator, 1. ) 
{
    check();
}

FloatVector::FloatVector( int dimension, const std::function<Float(int)>& generator, Float alpha )
: data( dimension )
{
    generatedatafrom( alpha, generator );
    check();
}



void FloatVector::check() const 
{
//     std::cout << "CHECK." << std::endl;
    assert( getdimension() >= 0 );
    if( getdimension() == 0 )
    std::cout << "WARNING: VECTOR OF DIMENSION ZERO." << std::endl;
    assert( getdimension() == data.size() );
}

void FloatVector::print( std::ostream& output ) const 
{
    check();
    output << "float vector of dimension: " << getdimension() << std::endl;
    for( int p = 0; p < getdimension(); p++ )
	output << p << ": " << getentry(p) << std::endl;
}

void FloatVector::printplain( std::ostream& output ) const 
{
    check();
    output << getdimension() << nl;
    for( int p = 0; p < getdimension(); p++ )
        output << getentry(p) << nl;
}



int FloatVector::getdimension() const 
{
    /* No check at this point */
    return data.size();
}

Float FloatVector::setentry( int p, Float value )
{
    check();
    assert( 0 <= p && p < data.size() );
    data.at(p) = value;
    return data.at(p);
}

Float FloatVector::getentry( int p ) const 
{
    check();
    assert( 0 <= p && p < data.size() );
    return data.at(p);
}

Float& FloatVector::at( int p )
{
    check();
    return (*this)[p];
}

const Float& FloatVector::at( int p ) const
{
    check();
    return (*this)[p];
}

Float& FloatVector::operator[]( int p )
{
    check();
    assert( 0 <= p && p < data.size() );
    return data.at(p);
}

const Float& FloatVector::operator[]( int p ) const
{
    check();
    assert( 0 <= p && p < data.size() );
    return data.at(p);
}
                
const std::vector<Float>& FloatVector::getdata() const
{
    check();
    return data;
}

        




void FloatVector::zero() 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
	setentry( p, 0. ); 
}

void FloatVector::random() 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, sqrt( rand() ) ); 
}

void FloatVector::scale( Float alpha ) 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, alpha * getentry( p ) ); 
}



FloatVector FloatVector::getslice( int base, int len ) const 
{
    check();
    assert( 0 <= base && base < getdimension() );
    assert( 0 <= len && len < getdimension() );
    assert( 0 <= base+len && base+len < getdimension() );
    FloatVector ret( len );
    for( int p = 0; p < len; p++ )
      ret[p] = data.at( base + p );
    return ret;
}

void FloatVector::setslice( int base, const FloatVector& slice )
{
    check();
    const int len = slice.getdimension();
    assert( 0 <= base && base < getdimension() );
    assert( 0 <= len && len < getdimension() );
    assert( 0 <= base+len && base+len < getdimension() );
    for( int p = 0; p < len; p++ )
      data.at( base + p ) = slice[p];
}
        
        








void FloatVector::copydatafrom( const FloatVector& source )
{
    check();
    copydatafrom( 1., source );
}

void FloatVector::copydatafrom( Float scaling, const FloatVector& source )
{
    check();
    assert( getdimension() == source.getdimension() );
        for( int p = 0; p < getdimension(); p++ )
            setentry( p, scaling * source.getentry( p ) );         
}
        
        
void FloatVector::generatedatafrom( const std::function<Float(int)>& generator )
{
    check();
    generatedatafrom( 1., generator );
}

void FloatVector::generatedatafrom( Float scaling, const std::function<Float(int)>& generator )
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, scaling * generator( p ) );         
}
        
        
void FloatVector::adddatafrom( const FloatVector& source )
{
    check();
    adddatafrom( 1., source );
}

void FloatVector::adddatafrom( Float scalingsource, const FloatVector& source )
{
    check();
    adddatafrom( 1., scalingsource, source );
}
        
void FloatVector::adddatafrom( Float scalingself, Float scalingsource, const FloatVector& source )
{
    check();
    assert( getdimension() == source.getdimension() );
        for( int p = 0; p < getdimension(); p++ )
            setentry( p, scalingself * getentry( p ) + scalingsource * source.getentry( p ) );         
}


Float FloatVector::scalarproductwith( const FloatVector& right ) const
{
    check();
    assert( getdimension() == right.getdimension() );
    Float ret = 0.;
    for( int p = 0; p < getdimension(); p++ )
        ret += getentry(p) * right.getentry(p);
    return ret;
}

Float FloatVector::norm() const 
{
    check();
    return sqrt( scalarproductwith( *this ) );
}

Float FloatVector::maxnorm() const
{
    check();
    assert( getdimension() > 0 );
    Float ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret = std::max( ret, absolute( data.at(d) ) );
    return ret;
}

Float FloatVector::lpnorm( Float p ) const
{
    check();
    assert( p > 0 );
    assert( std::isfinite(p) );
    assert( std::isnormal(p) );
    
    Float ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret += pow( absolute( data.at(d) ), p );
    return pow( ret, 1./p );
}





