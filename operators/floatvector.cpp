
#include <cassert>
#include <cmath>
#include <algorithm>
#include <new>
#include <ostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"






FloatVector::FloatVector( const FloatVector& src )
: dimension( src.getdimension() ), pointer( new (std::nothrow) Float[ src.getdimension() ] )
{
    assert( dimension >= 0 and pointer != nullptr );
    for( int p = 0; p < dimension; p++ ) pointer[p] = src.pointer[p];
    FloatVector::check();
}

FloatVector::FloatVector( FloatVector&& src )
: dimension( src.getdimension() ), pointer( src.pointer )
{
    assert( dimension >= 0 and pointer != nullptr and pointer == src.pointer and dimension == src.dimension );
    src.pointer = nullptr;
    FloatVector::check();
}

FloatVector& FloatVector::operator=( const FloatVector& vec )
{
    assert( pointer     != nullptr );
    assert( vec.pointer != nullptr );
    assert( getdimension() == vec.getdimension() );
    
    if( this != &vec ) {
        for( int p = 0; p < dimension; p++ ) pointer[p] = vec.pointer[p];
    }
    
    check();
    return *this;
}

FloatVector& FloatVector::operator=( FloatVector&& vec )
{
//     assert( pointer     != nullptr );
    assert( vec.pointer != nullptr );
    assert( getdimension() == vec.getdimension() );
    
    if( this != &vec ) {
        delete[] this->pointer;
        this->pointer = vec.pointer;
        vec.pointer = nullptr;
    }
    
    check();
    return *this;
}

FloatVector::~FloatVector()
{
    FloatVector::check();
    if( pointer != nullptr ){
        delete[] pointer;
    }
}













FloatVector::FloatVector( int dim, Float initialvalue )
: dimension(dim), pointer( new (std::nothrow) Float[dim] )
{
    assert( dimension >= 0 and pointer != nullptr );
    for( int p = 0; p < dimension; p++ ) pointer[p] = initialvalue;
    FloatVector::check();
}

FloatVector::FloatVector( const FloatVector& src, Float alpha )
: dimension( src.getdimension() ), pointer( new (std::nothrow) Float[ src.getdimension() ] )
{
    assert( dimension >= 0 and pointer != nullptr );
    for( int p = 0; p < dimension; p++ ) pointer[p] = alpha * src.pointer[p];
    FloatVector::check();
}

FloatVector::FloatVector( FloatVector&& src, Float alpha )
: dimension( src.getdimension() ), pointer( src.pointer )
{
    assert( dimension >= 0 and pointer != nullptr and pointer == src.pointer and dimension == src.dimension );
    src.pointer = nullptr;
    scale( alpha );
    FloatVector::check();
}

FloatVector::FloatVector( const std::vector<Float>& vals, Float alpha )
: FloatVector( SIZECAST( vals.size() ), [&vals](int i) -> Float{ return vals.at(i); }, alpha )
{
    FloatVector::check();
}

FloatVector::FloatVector( const std::vector<int>& vals, Float alpha )
: FloatVector( SIZECAST( vals.size() ), [&vals](int i) -> Float{ return vals.at(i); }, alpha )
{
    FloatVector::check();
}

FloatVector::FloatVector( const std::initializer_list<Float>& l )
: dimension( SIZECAST( l.size() ) ), pointer( new (std::nothrow) Float[l.size()] )
{
    assert( dimension >= 0 and pointer != nullptr );
    int i = 0;
    for( const Float& f : l ) pointer[i++] = f;
    FloatVector::check();
}

FloatVector::FloatVector( int dimension, const std::function<Float(int)>& generator, Float alpha )
: dimension( dimension ), pointer( new (std::nothrow) Float[dimension] )
{
    assert( dimension >= 0 and pointer != nullptr );
    generatedatafrom( alpha, generator );
    FloatVector::check();
}




















void FloatVector::check() const 
{
    #ifdef NDEBUG
    return;
    #endif
    assert( dimension >= 0 );
//     if( dimension > 0 ) assert( pointer != nullptr );
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

void FloatVector::lg() const 
{
    check();
    LOG << "float vector of dimension: " << getdimension() << std::endl;
    for( int p = 0; p < getdimension(); p++ )
        LOG << p << ": " << getentry(p) << std::endl;
}

void FloatVector::lgplain() const 
{
    check();
    LOG << getdimension() << std::endl;
    for( int p = 0; p < getdimension(); p++ )
        LOG << getentry(p) << nl;
}









FloatVector FloatVector::clone() const 
{
    return FloatVector( *this );
}











int FloatVector::getdimension() const 
{
    /* No check at this point */
    return dimension;
}

Float FloatVector::setentry( int p, Float value )
{
    //check();
    assert( 0 <= p && p < getdimension() );
    pointer[p] = value;
    return pointer[p];
}

Float FloatVector::getentry( int p ) const 
{
    //check();
    assert( 0 <= p && p < getdimension() );
    return pointer[p];
}












Float& FloatVector::at( int p ) &
{
    //check();
    assert( 0 <= p && p < getdimension() );
    return pointer[p];
}

const Float& FloatVector::at( int p ) const &
{
    //check();
    assert( 0 <= p && p < getdimension() );
    return pointer[p];
}

Float& FloatVector::operator[]( int p ) &
{
    //check();
    assert( 0 <= p && p < getdimension() );
    return pointer[p];
}

const Float& FloatVector::operator[]( int p ) const &
{
    //check();
    assert( 0 <= p && p < getdimension() );
    return pointer[p];
}
                
const std::vector<Float> FloatVector::getdata() const
{
    // check();
    std::vector<Float> ret( dimension );
    for( int p = 0; p < getdimension(); p++ ) 
        ret.at(p) = pointer[p];
       
    return ret;
}

        











void FloatVector::setentries( Float value ) 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, value ); 
}

void FloatVector::random() 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, gaussrand() ); 
}

void FloatVector::zero() 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, 0. ); 
}

void FloatVector::clear() 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, 0. ); 
}

void FloatVector::clear_if( const std::vector<bool>& mask ) 
{
    check();
    assert( mask.size() == getdimension() );
    for( int p = 0; p < getdimension(); p++ )
        if( mask[p] )
            setentry( p, 0. ); 
}

void FloatVector::clear_unless( const std::vector<bool>& mask ) 
{
    check();
    assert( mask.size() == getdimension() );
    for( int p = 0; p < getdimension(); p++ )
        if( not mask[p] )
            setentry( p, 0. ); 
}










FloatVector& FloatVector::normalize() 
{
    check();
    Float alpha = norm();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, getentry( p ) / alpha ); 
    return *this;
}

FloatVector& FloatVector::normalize( const LinearOperator& op ) 
{
    check();
    op.check();
    
    Float value = *this * ( op * (*this) );
    
    assert( std::isfinite(value) );
    assert( value >= 0 );
    
    value = sqrt(value);
    
    this->scaleinverse( value );
    
    return *this;
}

FloatVector& FloatVector::scale( Float alpha ) 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, alpha * getentry( p ) ); 
    return *this;
}

FloatVector& FloatVector::scaleinverse( Float alpha ) 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, getentry( p ) / alpha ); 
    return *this;
}

FloatVector& FloatVector::shift( Float delta ) 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, getentry( p ) + delta ); 
    return *this;
}

FloatVector& FloatVector::shiftnegative( Float delta ) 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, getentry( p ) - delta ); 
    return *this;
}





FloatVector FloatVector::getslice( int base, int len ) const 
{
    check();
    assert( 0 <= base && base < getdimension() );
    assert( 0 <= len && len <= getdimension() );
    assert( 0 <= base+len && base+len <= getdimension() );
    FloatVector ret( len );
    for( int p = 0; p < len; p++ )
      ret[p] = pointer[ base + p ];
    return ret;
}

void FloatVector::setslice( int base, int len, Float value ) 
{
    check();
    assert( 0 <= base && base < getdimension() );
    assert( 0 <= len && len <= getdimension() );
    assert( 0 <= base+len && base+len <= getdimension() );
    
    for( int p = 0; p < len; p++ ) pointer[ base + p ] = value;
}

void FloatVector::setslice( int base, const FloatVector& slice )
{
    check();
    const int len = slice.getdimension();
    assert( 0 <= base && base < getdimension() );
    assert( 0 <= len && len <= getdimension() );
    assert( 0 <= base+len && base+len <= getdimension() );
    for( int p = 0; p < len; p++ )
      pointer[ base + p ] = slice[p];
}

void FloatVector::addslice( int base, const FloatVector& slice, Float s )
{
    check();
    const int len = slice.getdimension();
    assert( 0 <= base && base < getdimension() );
    assert( 0 <= len && len <= getdimension() );
    assert( 0 <= base+len && base+len <= getdimension() );
    for( int p = 0; p < len; p++ )
      pointer[ base + p ] += s * slice[p];
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

Float FloatVector::scalarproductwith( const FloatVector& right, const std::vector<bool>& mask ) const
{
    check();
    assert( getdimension() == right.getdimension() );
    Float ret = 0.;
    for( int p = 0; p < getdimension(); p++ )
        if( not mask[p] )
            ret += getentry(p) * right.getentry(p);
    return ret;
}






Float FloatVector::average() const 
{
    return sum() / getdimension();
}

Float FloatVector::sum() const 
{
    check();
    assert( getdimension() > 0 );
    Float ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret = ret + pointer[d];
    return ret;
}

Float FloatVector::norm() const 
{
    return std::sqrt( norm_sq() );
}

Float FloatVector::norm_sq() const 
{
    check();
    Float ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret += absolute( pointer[d] ) * absolute( pointer[d] );
    return ret;
}

Float FloatVector::norm( const LinearOperator& op ) const 
{
    return std::sqrt( norm_sq( op ) );
}

Float FloatVector::norm_sq( const LinearOperator& op) const 
{
    check();
    op.check();
    Float ret = (*this) * ( op * (*this) );
    assert( std::isfinite(ret) );
    assert( ret >= 0 );
    return ret;
}

Float FloatVector::maxnorm() const
{
    check();
    assert( getdimension() > 0 );
    Float ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret = std::max( ret, absolute( pointer[d] ) );
    return ret;
}

Float FloatVector::sumnorm() const
{
    check();
    assert( getdimension() > 0 );
    Float ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret = ret + absolute( pointer[d] );
    return ret;
}

Float FloatVector::l2norm() const 
{
    return norm();
}

Float FloatVector::lpnorm( Float p ) const
{
    check();
    assert( p > 0 );
    assert( std::isfinite(p) );
    assert( std::isnormal(p) );
    
    Float ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret += power_numerical( absolute( pointer[d] ), p );
    return power_numerical( ret, 1./p );
}








bool FloatVector::isfinite() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( not std::isfinite( pointer[d] ) )
            return false;
    return true;
}

bool FloatVector::iszero() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( pointer[d] != 0. )
            return false;
    return true;
}

bool FloatVector::ispositive() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( not ( pointer[d] > 0. ) )
            return false;
    return true;
}

bool FloatVector::isnegative() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( not ( pointer[d] < 0. ) )
            return false;
    return true;
}

bool FloatVector::isnonnegative() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( not ( pointer[d] >= 0. ) )
            return false;
    return true;
}

bool FloatVector::isnonpositive() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( not ( pointer[d] <= 0. ) )
            return false;
    return true;
}



bool FloatVector::issmall( Float eps ) const 
{
    check();
    return this->norm() < eps;
}





Float* FloatVector::raw()
{
    return pointer;
}

const Float* FloatVector::raw() const
{
    return pointer;
}










