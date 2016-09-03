#ifndef INCLUDEGUARD_LINEAROPERATOR
#define INCLUDEGUARD_LINEAROPERATOR




#include <ostream>

#include "../basic.hpp"
#include "floatvector.hpp"



/******************
*** 
***  Linear Operator on an abstract level 
***  - applyAdd is a purely virtual function
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
        void apply( FloatVector& dest, const FloatVector& src, Float scaling = 1. ) const;
        
        /* x := s x + t A y */
        virtual void applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const = 0;
        
    private:
        
        int dimout;
        int dimin;
    
};
  
  
  
inline FloatVector operator*( const LinearOperator& op, const FloatVector& vec )
{
    FloatVector ret( op.getdimout() );
    op.apply( ret, vec );
    return ret;
}
  
inline std::ostream& operator<<( std::ostream& os, const LinearOperator& op )
{
    op.print( os );
    return os;
}
  

  
  
  
  
  
#endif