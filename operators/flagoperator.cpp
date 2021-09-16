
#include "flagoperator.hpp"

#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"







FlagOperator::FlagOperator( const LinearOperator& op, const std::vector<bool>& destflag, const std::vector<bool>& srcflag )
: LinearOperator( op.getdimout() , op.getdimin() ), 
  op( op ),
  destflag( destflag ),
  srcflag( srcflag )
{
    FlagOperator::check();
}

FlagOperator::FlagOperator( const LinearOperator& op, const std::vector<bool>& flag )
: FlagOperator( op, flag, flag )
{
    FlagOperator::check();
}

FlagOperator::~FlagOperator()
{
        /* Nothing */ 
}

void FlagOperator::check() const  
{
    LinearOperator::check();    
    assert( getdimin()  == srcflag.size()  );
    assert( getdimout() == destflag.size() );
}

std::string FlagOperator::text() const  
{
    return "Flag Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin());
}

void FlagOperator::print( std::ostream& os ) const  
{
    os << text() << std::endl;
}



std::vector<bool>& FlagOperator::getsrcflag()
{
    return srcflag;
}

const std::vector<bool>& FlagOperator::getsrcflag() const
{
    return srcflag;
}

std::vector<bool>& FlagOperator::getdestflag()
{
    return destflag;
}

const std::vector<bool>& FlagOperator::getdestflag() const
{
    return destflag;
}




void FlagOperator::apply( FloatVector& dest, const FloatVector& src, Float scaling ) const 
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    FloatVector cloned_src = src;
    cloned_src.clear_if( srcflag );
    
    op.apply( dest, cloned_src, scaling );
    
    dest.clear_if( destflag );
    
}
