
#include "sumoperator.hpp"

#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


SumOperator::SumOperator( const LinearOperator& left, const LinearOperator& right )
: SumOperator( left, right, 1., 1. )
{
    check();
}

SumOperator::SumOperator( const LinearOperator& left, const LinearOperator& right,
 Float leftscale, Float rightscale )
: LinearOperator(left.getdimout(),right.getdimin()), 
  left(left), right(right),
  leftscale(leftscale), rightscale(rightscale)
{
    check();
}
		
SumOperator::~SumOperator()
{
	/* Nothing */ 
}

void SumOperator::check() const  
{
    LinearOperator::check();
    left.check();
    right.check();
    assert( left.getdimin() == right.getdimin() );
    assert( left.getdimout() == right.getdimout() );
}

void SumOperator::print( std::ostream& os ) const  
{
    os << "Print Sum Operator: " << std::endl;
}



void SumOperator::applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const 
{
    check();
    dest.check();
    add.check();
    
    assert( getdimin() == getdimout() );
    assert( getdimout() == dest.getdimension() );
    assert( getdimin() == add.getdimension() );
      
    dest = s * dest + t * leftscale * ( left * add ) + t * rightscale * ( right * add );
    
    dest.check();
}