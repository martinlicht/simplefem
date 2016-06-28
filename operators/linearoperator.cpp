

#include "linearoperator.hpp"

    LinearOperator::LinearOperator( int out, int in )
    : dimout( out ), dimin( in )
    {
		/* Nothing left to do */
    }
    
    LinearOperator::~LinearOperator()
    {
    
      /* todo */
      
    }
    
    int LinearOperator::getdimout() const
    {
      return dimout;
    }
    
    int LinearOperator::getdimin() const
    {
      return dimin;
    }
    
    void LinearOperator::check() const
    {
	  assert( dimout >= 0 && dimin >= 0 );
    }
    
    void LinearOperator::print( std::ostream& os ) const
    {
      os << "Linear operator" << std::endl;
    }
    
    /* x := A y */
    void LinearOperator::apply( FloatVector& dest, const FloatVector& src, Float f ) const
    {
      
      applyadd( dest, src, 0.0, f );
      
    }
    
    
    
  
  
