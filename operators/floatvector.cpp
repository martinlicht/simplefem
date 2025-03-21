
#include <cmath>
#include <cstddef>

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <new>
#include <string>
#include <vector>

#include "../base/include.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"
#include "../utility/random.hpp"






FloatVector::FloatVector( const FloatVector& src )
: dimension( src.dimension ), pointer( new (std::nothrow) Float[ src.dimension ] )
{
    assert( dimension >= 0 and pointer != nullptr );
    for( int p = 0; p < dimension; p++ ) pointer[p] = src.pointer[p];
    FloatVector::check();
}

FloatVector::FloatVector( FloatVector&& src ) noexcept
: dimension( std::move(src.dimension) ), pointer( std::move(src.pointer) )
{
    assert( dimension >= 0 and pointer != nullptr and pointer == src.pointer and dimension == src.dimension );
    src.pointer = nullptr;
    FloatVector::check();
}

FloatVector& FloatVector::operator=( const FloatVector& vec )
{
    assert( dimension == vec.dimension );
    
    if( this != &vec ) {
        assert( vec.pointer != nullptr );
        for( int p = 0; p < dimension; p++ ) pointer[p] = vec.pointer[p];
    }
    
    check();
    return *this;
}

FloatVector& FloatVector::operator=( FloatVector&& vec ) noexcept
{
    assert( dimension == vec.dimension );
    
    if( this != &vec ) {
        // std::swap( this->pointer, vec.pointer );
        delete[] this->pointer;
        this->pointer = vec.pointer;
        vec.pointer = nullptr;
    }
    
    check();
    return *this;
}

FloatVector::~FloatVector() noexcept
{
    if( pointer != nullptr ){
        FloatVector::check();
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

FloatVector::FloatVector( const FloatVector& src, Float scaling )
: dimension( src.dimension ), pointer( new (std::nothrow) Float[ src.dimension ] )
{
    assert( dimension >= 0 and pointer != nullptr );
    for( int p = 0; p < dimension; p++ ) pointer[p] = scaling * src.pointer[p];
    FloatVector::check();
}

FloatVector::FloatVector( FloatVector&& src, Float scaling )
: dimension( src.dimension ), pointer( src.pointer )
{
    assert( dimension >= 0 and pointer != nullptr and pointer == src.pointer and dimension == src.dimension );
    src.pointer = nullptr;
    scale( scaling );
    FloatVector::check();
}

FloatVector::FloatVector( const std::vector<Float>& vals, Float scaling )
: FloatVector( SIZECAST( vals.size() ), [&vals](int i) -> Float{ return vals.at(i); }, scaling )
{
    FloatVector::check();
}

FloatVector::FloatVector( const std::vector<int>& vals, Float scaling )
: FloatVector( SIZECAST( vals.size() ), [&vals](int i) -> Float{ return vals.at(i); }, scaling )
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

FloatVector::FloatVector( int dimension, const std::function<Float(int)>& generator, Float scaling )
: dimension( dimension ), pointer( new (std::nothrow) Float[dimension] )
{
    assert( dimension >= 0 and pointer != nullptr );
    generatedatafrom( scaling, generator );
    FloatVector::check();
}




















void FloatVector::check() const 
{
    #ifdef NDEBUG
    return;
    #endif
    assert( dimension >= 0 );
    // if( pointer == nullptr ) assert( dimension == 0 );
    // if( dimension > 0 ) assert( pointer != nullptr );
}

std::string FloatVector::text() const 
{
    check();
    std::string ret = "float vector of dimension: " + std::to_string( getdimension() );
    for( int p = 0; p < getdimension(); p++ ) {
        // ret = ret + "\n" + std::to_string(p) + ": " + std::to_string(getentry(p));
        ret = ret + "\n" + std::to_string(p) + ": " + printf_into_string( "% .17le", (double)(safedouble)getentry(p) );
    }
    return ret;
}

std::string FloatVector::data_as_text( bool indexed, bool print_rowwise ) const
{
    const int nc_precision = 10;

    const int nc_width = 7 + nc_precision;
    
    std::string ret;

    if( print_rowwise ){
        if(indexed) for( int p = 0; p < getdimension(); p++ ) ret += printf_into_string(    "%*d ", nc_width, p );
        if(indexed) ret += '\n';
        for( int p = 0; p < getdimension(); p++ ) ret += printf_into_string( "%*.*le ", nc_width, nc_precision, (double)(safedouble)getentry(p) );
    } else {
        for( int p = 0; p < getdimension(); p++ )
            if(indexed) 
                ret += printf_into_string( "%*d : %*.*le\n", nc_width, p, nc_width, nc_precision, (double)(safedouble)getentry(p) );
            else 
                ret += printf_into_string( "%*.*le\n", nc_width, nc_precision, (double)(safedouble)getentry(p) );
    }

    return ret;
}


// void FloatVector::print( std::ostream& output ) const 
// {
//     check();
//     output << "float vector of dimension: " << getdimension() << nl;
//     for( int p = 0; p < getdimension(); p++ )
//         output << p << ": " << getentry(p) << nl;
// }









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
                
std::vector<Float> FloatVector::getdata() const
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
        setentry( p, gaussian_variable() ); 
}

void FloatVector::random_within_range( Float min, Float max )
{
    check();
    for( int p = 0; p < getdimension(); p++ ) {
        Float value = (max-min)*random_uniform() + min; 
        assert( min <= value and value <= max );
        setentry( p, value ); 
    }   
}
        
void FloatVector::to_absolute() 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, absolute( getentry(p) ) ); 
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
    
    Float value_sqrt = std::sqrt(value);
    
    this->scaleinverse( value_sqrt );
    
    return *this;
}

FloatVector& FloatVector::scale( Float factor ) 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, factor * getentry( p ) ); 
    return *this;
}

FloatVector& FloatVector::scaleinverse( Float divisor ) 
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, getentry( p ) / divisor ); 
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

FloatVector& FloatVector::reciprocal()
{
    check();
    for( int p = 0; p < getdimension(); p++ )
        setentry( p, 1./getentry( p ) ); 
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

void FloatVector::setslice( int base, int len, Float uniform_value ) 
{
    check();
    assert( 0 <= base && base < getdimension() );
    assert( 0 <= len && len <= getdimension() );
    assert( 0 <= base+len && base+len <= getdimension() );
    
    for( int p = 0; p < len; p++ ) pointer[ base + p ] = uniform_value;
}

void FloatVector::setslice( int base, const FloatVector& source )
{
    check();
    const int len = source.getdimension();
    assert( 0 <= base && base < getdimension() );
    assert( 0 <= len && len <= getdimension() );
    assert( 0 <= base+len && base+len <= getdimension() );
    for( int p = 0; p < len; p++ )
        pointer[ base + p ] = source[p];
}

void FloatVector::addslice( int base, const FloatVector& summand, Float scaling )
{
    check();
    const int len = summand.getdimension();
    assert( 0 <= base && base < getdimension() );
    assert( 0 <= len && len <= getdimension() );
    assert( 0 <= base+len && base+len <= getdimension() );
    for( int p = 0; p < len; p++ )
        pointer[ base + p ] += scaling * summand[p];
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
        
        
void FloatVector::adddatafrom( const FloatVector& summand )
{
    check();
    adddatafrom( 1., summand );
}

void FloatVector::adddatafrom( Float scalingsummand, const FloatVector& summand )
{
    check();
    adddatafrom( 1., scalingsummand, summand );
}
        
void FloatVector::adddatafrom( Float scalingdest, Float scalingsummand, const FloatVector& summand )
{
    check();
    assert( getdimension() == summand.getdimension() );
        for( int p = 0; p < getdimension(); p++ )
            setentry( p, scalingdest * getentry( p ) + scalingsummand * summand.getentry( p ) );         
}


Float FloatVector::scalarproductwith( const FloatVector& right ) const
{
    check();
    assert( getdimension() == right.getdimension() );
    long double ret = 0.;
    for( int p = 0; p < getdimension(); p++ )
        ret += getentry(p) * right.getentry(p);
    return static_cast<Float>(ret);
}

Float FloatVector::scalarproductwith( const FloatVector& right, const std::vector<bool>& mask ) const
{
    check();
    assert( getdimension() == right.getdimension() );
    long double ret = 0.;
    for( int p = 0; p < getdimension(); p++ )
        if( not mask[p] )
            ret += getentry(p) * right.getentry(p);
    return static_cast<Float>(ret);
}






Float FloatVector::min() const 
{
    check();
    assert( getdimension() > 0 );
    Float ret = pointer[0];
    for( int d = 1; d < getdimension(); d++ )
        ret = ( ret < pointer[d] ) ? ret : pointer[d];
    return ret;
}

Float FloatVector::average() const 
{
    assert( getdimension() > 0 );
    return sum() / getdimension();
}

Float FloatVector::sum() const 
{
    check();
    long double ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret = ret + pointer[d];
    return static_cast<Float>(ret);
}

Float FloatVector::norm() const 
{
    return std::sqrt( norm_sq() );
}

Float FloatVector::norm_sq() const 
{
    check();
    long double ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret += absolute( pointer[d] * pointer[d] );
    return static_cast<Float>(ret);
}

Float FloatVector::norm( const LinearOperator& op ) const 
{
    return std::sqrt( norm_sq( op ) );
}

Float FloatVector::norm_sq( const LinearOperator& op) const 
{
    check();
    op.check();
    long double ret = (*this) * ( op * (*this) );
    assert( std::isfinite(ret) );
    assert( ret >= 0 );
    return static_cast<Float>(ret);
}

Float FloatVector::maxnorm() const
{
    check();
    assert( getdimension() > 0 );
    Float ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret = maximum( ret, absolute( pointer[d] ) );
    return ret;
}

Float FloatVector::sumnorm() const
{
    check();
    assert( getdimension() > 0 );
    long double ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret = ret + absolute( pointer[d] );
    return static_cast<Float>(ret);
}

Float FloatVector::l2norm() const 
{
    return norm();
}

Float FloatVector::lpnorm( Float p, Float inner_weight ) const
{
    check();
    assert( p > 0 );
    assert( std::isfinite(p) );
    assert( std::isnormal(p) );
    assert( std::isfinite(inner_weight) and inner_weight > 0. );
    
    long double ret = 0.;
    for( int d = 0; d < getdimension(); d++ )
        ret += power_numerical( absolute( pointer[d] ), p );
    ret *= inner_weight;
    return power_numerical( static_cast<Float>(ret), 1./p );
}








bool FloatVector::is_finite() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( not std::isfinite( pointer[d] ) )
            return false;
    return true;
}

bool FloatVector::is_zero() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( pointer[d] != 0. )
            return false;
    return true;
}

bool FloatVector::is_positive() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( not ( pointer[d] > 0. ) )
            return false;
    return true;
}

bool FloatVector::is_negative() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( not ( pointer[d] < 0. ) )
            return false;
    return true;
}

bool FloatVector::is_nonnegative() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( not ( pointer[d] >= 0. ) )
            return false;
    return true;
}

bool FloatVector::is_nonpositive() const 
{
    check();
    for( int d = 0; d < getdimension(); d++ )
        if( not ( pointer[d] <= 0. ) )
            return false;
    return true;
}



bool FloatVector::is_numerically_small( Float threshold ) const 
{
    check();
    return this->norm() < threshold;
}





Float* FloatVector::raw()
{
    return pointer;
}

const Float* FloatVector::raw() const
{
    return pointer;
}



/* Memory size */
        
std::size_t FloatVector::memorysize() const
{
    return sizeof(*this) + this->dimension * sizeof(decltype(*pointer));
}


bool FloatVector::is_equal_to( const FloatVector& vector_left, const FloatVector& vector_right )
{
    assert( vector_left.getdimension() == vector_right.getdimension() );
    for( int i = 0; i < vector_left.getdimension(); i++ )
        if( vector_left[i] != vector_right[i] )
            return false;
    return true;
}

        




FloatVector interlace( const FloatVector& first, const FloatVector& second )
{
    Assert( first.getdimension() == second.getdimension() );
    FloatVector ret( first.getdimension() * 2 );
    for( int i = 0; i < first.getdimension(); i++ )
    {
        ret[2*i+0] = first[i];
        ret[2*i+1] = second[i];
    }
    return ret;
}
