
#include "sumoperator.hpp"

#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


SumOperator::SumOperator( const LinearOperator& left, const LinearOperator& right,
 Float leftscale, Float rightscale )
: LinearOperator(left.getdimout(),right.getdimin()), 
  left(left), right(right),
  leftscale(leftscale), rightscale(rightscale)
{
    SumOperator::check();
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



void SumOperator::apply( FloatVector& dest, const FloatVector& src, Float scaling ) const 
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    dest = scaling * leftscale * ( left * src ) + scaling * rightscale * ( right * src );
}



