
#include "diagonaloperator.hpp"

#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


DiagonalOperator::DiagonalOperator( int dimension, Float scale )
: LinearOperator(dimension,dimension), 
  diagonal( FloatVector( dimension, scale ) )
{
    DiagonalOperator::check();
}

DiagonalOperator::DiagonalOperator( int dimension, const FloatVector& dia )
: LinearOperator(dimension,dimension), diagonal(dia)
{
    DiagonalOperator::check();
}

DiagonalOperator::DiagonalOperator( int dimension, const ScalingOperator& scaling )
: LinearOperator(dimension,dimension), 
  diagonal( FloatVector( dimension, scaling.getscaling() ) )
{
    DiagonalOperator::check();
}

DiagonalOperator::~DiagonalOperator()
{
        /* Nothing */ 
}

void DiagonalOperator::check() const  
{
    LinearOperator::check();    
    diagonal.check();
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


void DiagonalOperator::apply( FloatVector& dest, const FloatVector& src, Float scaling ) const 
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == getdimout() );
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    for( int p = 0; p < getdimout(); p++ ) {
        dest.setentry( p, scaling * diagonal.at(p) * src.getentry( p ) );
    }
    
}
