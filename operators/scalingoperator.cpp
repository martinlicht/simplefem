
#include <iostream>
#include "../basic.hpp"
#include "scalingoperator.hpp"


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

Float ScalingOperator::getscaling()
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
		


void ScalingOperator::applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const 
{
	
	assert( getdimin() == getdimout() );
	assert( getdimout() == dest.getdimension() );
	assert( getdimin() == add.getdimension() );
	
	for( int p = 0; p < getdimin(); p++ ) {
		dest.setentry( p, 
			s * dest.getentry( p ) + t * scaling * add.getentry( p )
		);
	}
	
	
}