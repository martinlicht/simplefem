
#include "blockoperator.hpp"

#include <iostream>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


BlockOperator::BlockOperator( int dimension_out, int dimension_in, const std::vector<std::vector<LinearOperator*>>& ops )
: LinearOperator( dimension_out, dimension_in ), ops(ops)
{
  
}

BlockOperator::~BlockOperator()
{
    /* Nothing */ 
}

void BlockOperator::check() const  
{
  std::clog << "Check Block Operator" << std::endl;
  LinearOperator::check();
  
  if( ops.size() == 0 ) {
    
    std::clog << "zero rows in block operator" << std::endl;
    
  } else if ( ops.at(0).size() == 0 ) {
    
    std::clog << "empty rows in block operator" << std::endl;
    
    for( int r = 1; r < ops.size(); r++ )
      assert( ops.at(r).size() == ops.at(0).size() ); 
    
  } else {
    
    std::clog << "full operator block matrix" << std::endl;
    
    int temp_out = 0;
    
    for( int r = 0; r < ops.size(); r++ ) {
      
      assert( ops.size() > 0 );
      assert( ops.at(r).size() > 0 );
      
      assert( ops.at(r).size() == ops.at(0).size() );
      for( int c = 1; c < ops.at(0).size(); c++ )
        assert( ops.at(r).at(0)->getdimout() == ops.at(r).at(c)->getdimout() ); 
      
      int temp_in = 0;
      for( int c = 0; c < ops.at(0).size(); c++ )
        temp_in += ops.at(r).at(c)->getdimout(); 
      assert( temp_in == getdimin() );
        
      temp_out += ops.at(r).at(0)->getdimout();
      
    }
    
    assert( temp_out == getdimout() );
    
    std::clog << "full operator block matrix" << std::endl;
    
    for( int c = 0; c < ops.at(0).size(); c++ ) {
      
      for( int r = 0; r < ops.size(); r++ ) {
        
        assert( ops.at(r).at(c)->getdimin() == ops.at(0).at(c)->getdimin() );
        
      }
      
    }
    
  }
  
}

void BlockOperator::print( std::ostream& os ) const  
{
    os << "Print Block Operator" << std::endl;
}





FloatVector BlockOperator::apply( const FloatVector& src, Float s ) const 
{
    
    check();
    src.check();
    assert( getdimin() == src.getdimension() );
    
    FloatVector ret( getdimout() );
    
    int base_out = 0;
    
    for( int r = 0; r < ops.size(); r++ ) {
      
      int base_in  = 0;
    
      for( int c = 0; c < ops.at(0).size(); c++ ) {
        
        ret.addslice(
          base_out,
          ops.at(r).at(c)->apply( src.getslice( base_in, ops.at(r).at(c)->getdimin() ), s ),
          1.
        );
        
        base_in += ops.at(r).at(c)->getdimin();
        
      }
      
      assert( base_in == getdimin() );
      
      base_out += ops.at(r).at(0)->getdimout();
      
    }
    
    assert( base_out == getdimout() );
    
    return ret;
   
}


