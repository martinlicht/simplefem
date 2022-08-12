
#include "complexoperator.hpp"

#include <cmath>
#include <functional>
#include <ostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"
#include "simpleoperators.hpp"


FloatVector RealPart( const FloatVector& vec )
{
    assert( vec.getdimension() % 2 == 0 );
    return vec.getslice( 0, vec.getdimension() / 2 );
}

FloatVector ImaginaryPart( const FloatVector& vec )
{
    assert( vec.getdimension() % 2 == 0 );
    return vec.getslice( vec.getdimension() / 2, vec.getdimension() / 2 );
}

FloatVector ComplexFloatVector( const FloatVector& real, const FloatVector& imag )
{
    assert( real.getdimension() == imag.getdimension() );
    FloatVector ret( 2 * real.getdimension() );
    ret.setslice(                   0, real );
    ret.setslice( real.getdimension(), imag );
    return ret;
}





ComplexOperator::ComplexOperator( const LinearOperator& real, const LinearOperator& imag )
: LinearOperator( real.getdimout() + imag.getdimout(), real.getdimin() + imag.getdimin() ), 
part_real( real ), part_imag( imag ) 
{
    ComplexOperator::check();
}

ComplexOperator::~ComplexOperator()
{
    /* Nothing */ 
}

void ComplexOperator::check() const  
{
    LinearOperator::check();    
    assert( part_real.getdimin()  == part_imag.getdimin()  );
    assert( part_real.getdimout() == part_imag.getdimout() );
}

std::string ComplexOperator::text() const 
{
    return text( false ); // TODO use embellish...
}

std::string ComplexOperator::text( const bool embellish ) const 
{
    return "Complex operator";
}


const LinearOperator& ComplexOperator::real() const 
{
    return part_real;
}

const LinearOperator& ComplexOperator::imaginary() const 
{
    return part_imag;
}

void ComplexOperator::apply( FloatVector& dest, const FloatVector& src, Float scaling ) const 
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == getdimout() );
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    auto src_real = RealPart( src );
    auto src_imag = ImaginaryPart( src );

    auto dest_real = scaling * ( part_real * src_real - part_imag * src_imag );
    auto dest_imag = scaling * ( part_real * src_imag + part_imag * src_real );

    dest = ComplexFloatVector( dest_real, dest_imag );
}

