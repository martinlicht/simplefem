
#include "blockdiagonaloperator.hpp"

#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


BlockDiagonalOperator::BlockDiagonalOperator( int dimension_out, int dimension_in, const std::vector<LinearOperator*>& ops )
: LinearOperator( dimension_out, dimension_in ), ops(ops)
{}

BlockDiagonalOperator::~BlockDiagonalOperator()
{
    /* Nothing */ 
}

void BlockDiagonalOperator::check() const  
{
    LinearOperator::check();
    
    int temp_in = 0;
    int temp_out = 0;
    
    for( const auto& op : ops ) { 
      op->check();
      temp_in += op->getdimin();
      temp_out += op->getdimout();
    }
    
    assert( temp_in == getdimin() );
    assert( temp_out == getdimout() );
}

void BlockDiagonalOperator::print( std::ostream& os ) const  
{
    os << "Print Block Diagonal Operator" << std::endl;
}





void BlockDiagonalOperator::apply( FloatVector& dest, const FloatVector& src, Float scaling ) const 
{
    check();
    src.check();
    dest.check();

    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );

    int base_in = 0;
    int base_out = 0;
    
    for( const auto& op : ops ) {
        
        dest.setslice( base_out,
                      op->apply( src.getslice( base_in, op->getdimin() ), scaling )
                    );
        
        base_in += op->getdimin();
        base_out += op->getdimout();
    }
    
    assert( base_in == getdimin() );
    assert( base_out == getdimout() );
    
}
