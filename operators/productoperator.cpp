
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



void ProductOperator::applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const 
{
    check();
    dest.check();
    add.check();
    
    assert( getdimin() == getdimout() );
    assert( getdimout() == dest.getdimension() );
    assert( getdimin() == add.getdimension() );
      
    FloatVector temp = right * add;
    dest = s * dest + t * ( left * temp );
    
    dest.check();
}