#ifndef INCLUDEGUARD_OPERATOR_ABSTRACT
#define INCLUDEGUARD_OPERATOR_ABSTRACT




#include <ostream>

#include "../basic.hpp"
#include "floatvector.hpp"



/******************
*** 
***  Linear Operator on an abstract level 
***  - applyAdd is a purely virtual function that implements y <- A x with optional scaling of the result
***  - all other apply adds are defined in terms of that, but overwriting in derived classes is possible
***
******************/

class LinearOperator
{

    public:
        
        explicit LinearOperator( int, int );
        virtual ~LinearOperator();
        
        int getdimin() const;
        int getdimout() const;
        
        virtual void check() const;
        virtual void print( std::ostream& os ) const;
        
        /* Apply the operator */
        
        /* x := s A y */
        virtual FloatVector apply( const FloatVector& src, Float scaling = 1. ) const = 0;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling = 1. ) const;
        
        /* x := s x + t A y */
        virtual void applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const;
        
    private:
        
        int dimout;
        int dimin;
    
};
  
  
  
inline FloatVector operator*( const LinearOperator& op, const FloatVector& vec )
{
    op.check();
    vec.check();
    FloatVector ret( op.getdimout() );
    op.apply( ret, vec );
    vec.check();
    return ret;
}
  
inline std::ostream& operator<<( std::ostream& os, const LinearOperator& op )
{
    op.check();
    op.print( os );
    op.check();
    return os;
}
  

  
  
  
  
  
#endif