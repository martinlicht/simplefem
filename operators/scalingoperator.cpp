
#include "scalingoperator.hpp"

#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


ScalingOperator::ScalingOperator( int dimension, Float s )
: LinearOperator(dimension,dimension), scaling(s)
{}

ScalingOperator::~ScalingOperator()
{
    /* Nothing */ 
}

void ScalingOperator::check() const  
{
    LinearOperator::check();
}

void ScalingOperator::print( std::ostream& os ) const  
{
    os << "Print Scaling Operator with scaling: " << scaling << std::endl;
}


void ScalingOperator::setscaling( Float s )
{
    scaling = s;
}

Float ScalingOperator::getscaling() const
{
    return scaling;
}





void ScalingOperator::apply( FloatVector& dest, const FloatVector& src, Float s ) const 
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == getdimout() );
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    for( int p = 0; p < getdimin(); p++ )
        dest.setentry( p, s * scaling * src.getentry( p ) );
        
}
