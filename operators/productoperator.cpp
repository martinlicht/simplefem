
#include "productoperator.hpp"

#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


ProductOperator::ProductOperator( const LinearOperator& left, const LinearOperator& right )
: LinearOperator(left.getdimout(),right.getdimin()), 
  left(left), right(right)
{
    check();
}

ProductOperator::~ProductOperator()
{
    /* Nothing */ 
}

void ProductOperator::check() const  
{
    LinearOperator::check();
    left.check();
    right.check();
    assert( left.getdimin() == right.getdimout() );
    assert( getdimin() == right.getdimin() );
    assert( getdimout() == left.getdimout() );
}

void ProductOperator::print( std::ostream& os ) const  
{
    os << "Print Product Operator: " << std::endl;
}



FloatVector ProductOperator::apply( const FloatVector& src, Float scaling ) const 
{
    check();
    src.check();
    assert( getdimin() == src.getdimension() );
    
    return scaling * ( left * ( right * src ) );
}







