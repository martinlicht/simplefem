
#include "simpleoperators.hpp"

#include <cmath>
#include <functional>
#include <ostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


ScalingOperator::ScalingOperator( int dimension, Float s )
: LinearOperator( dimension ), 
  scaling(s)
{
    ScalingOperator::check();
}

ScalingOperator::~ScalingOperator()
{
    ScalingOperator::check();
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














DiagonalOperator::DiagonalOperator( int dimension, Float scale )
: LinearOperator( dimension ), 
  diagonal( FloatVector( dimension, scale ) )
{
    DiagonalOperator::check();
}

DiagonalOperator::DiagonalOperator( const FloatVector& dia )
: LinearOperator( dia.getdimension() ), 
  diagonal( dia )
{
    DiagonalOperator::check();
}

DiagonalOperator::DiagonalOperator( FloatVector&& dia )
: LinearOperator( dia.getdimension() ), 
  diagonal( std::move(dia) )
{
    DiagonalOperator::check();
}

DiagonalOperator::DiagonalOperator( int dimension, const ScalingOperator& scaling )
: LinearOperator( dimension ), 
  diagonal( FloatVector( dimension, scaling.getscaling() ) )
{
    DiagonalOperator::check();
}

DiagonalOperator::DiagonalOperator( int dimension, const std::function<Float(int)>& generator )
: LinearOperator( dimension ), 
  diagonal( FloatVector( dimension, generator ) )
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

const DiagonalOperator DiagonalOperator::sqrt() const
{
    return DiagonalOperator( getdimin(), [this](int i) -> Float {
        assert( this->diagonal[i] >= 0 ); 
        return std::sqrt( diagonal[i] );
    } );
}

