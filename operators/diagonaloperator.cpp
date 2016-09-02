
#include "diagonaloperator.hpp"

#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


DiagonalOperator::DiagonalOperator( int dimension, const FloatVector& dia )
: LinearOperator(dimension,dimension), diagonal(dia)
{
    check();
}
		
DiagonalOperator::DiagonalOperator( int dimension, const ScalingOperator& scaling )
: LinearOperator(dimension,dimension), 
  diagonal( FloatVector( dimension, scaling.getscaling() ) )
{
    check();
}
		
DiagonalOperator::~DiagonalOperator()
{
	/* Nothing */ 
}

void DiagonalOperator::check() const  
{
	std::cout << "Check Diagonal Operator:" << std::endl;
	diagonal.check();
	LinearOperator::check();	
    assert( getdimin() == getdimout() );
    assert( getdimin() == diagonal.getdimension() );
}

void DiagonalOperator::print( std::ostream& os ) const  
{
	os << "Print Diagonal Operator with diagonal: " 
	   << diagonal << std::endl;
}


FloatVector& DiagonalOperator::getdiagonal()
{
	return diagonal;
}

const FloatVector& DiagonalOperator::getdiagonal() const
{
	return diagonal;
}


void DiagonalOperator::applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const 
{
	assert( getdimin() == getdimout() );
	assert( getdimout() == dest.getdimension() );
	assert( getdimin() == add.getdimension() );
	
	for( int p = 0; p < getdimin(); p++ ) {
		dest.setentry( p, 
			s * dest.getentry( p ) + t * diagonal.getentry( p ) * add.getentry( p )
		);
	}	
}