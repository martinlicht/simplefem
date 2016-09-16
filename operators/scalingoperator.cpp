
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
    std::cout << "Check Scaling Operator with scaling: " << scaling << std::endl;
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



// SparseMatrix ScalingOperator::cloneSparseMatrix() const 
// {
	// SparseMatrix M(dimin,dimout);
	// for( int p = 0; p < dimin; p++ )
		// M.addentry( p, p, scaling );
	// return M;
// };
		


FloatVector ScalingOperator::apply( const FloatVector& src, Float scaling ) const 
{
    check();
    src.check();
    
    assert( getdimin() == getdimout() );
    assert( getdimin() == src.getdimension() );
    
    FloatVector ret( getdimout() );
    
    for( int p = 0; p < getdimin(); p++ )
        ret.setentry( p, scaling * src.getentry( p ) );
    
    return ret;
    
}