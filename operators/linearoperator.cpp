

#include "linearoperator.hpp"

#include "floatvector.hpp"

    LinearOperator::LinearOperator( int out, int in )
    : dimout( out ), dimin( in )
    {
      check();
    }
    
    LinearOperator::~LinearOperator()
    {
    
      check();
      
    }
    
    int LinearOperator::getdimout() const
    {
      /* No check here */
      return dimout;
    }
    
    int LinearOperator::getdimin() const
    {
      /* No check here */
      return dimin;
    }
    
    void LinearOperator::check() const
    {
      assert( dimout >= 0 && dimin >= 0 );
    }
    
    void LinearOperator::print( std::ostream& os ) const
    {
      check();
      os << "Linear operator" << std::endl;
    }
    
    /* x := A y */
    void LinearOperator::apply( FloatVector& dest, const FloatVector& src, Float scaling ) const
    {
      check();
      dest.copydatafrom( apply( src, scaling ) );
    }
    
    /* x := s x + t A y */
    void LinearOperator::applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const
    {
      FloatVector temp = apply( add, t );
      dest.adddatafrom( s, 1., temp );
    }
        
    
    
    
  
  
