#ifndef INCLUDEGUARD_OPERATOR_LINEAROPERATOR_HPP
#define INCLUDEGUARD_OPERATOR_LINEAROPERATOR_HPP




#include <memory>
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
        
        explicit LinearOperator() = delete;
        
        explicit LinearOperator( int, int );

        
        explicit LinearOperator( const LinearOperator& )       = default;
        explicit LinearOperator( LinearOperator&& )            = default;
        LinearOperator& operator=( const LinearOperator& vec ) = default;
        LinearOperator& operator=( LinearOperator&& vec )      = default;
        
        virtual ~LinearOperator();

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const&
        {
            unreachable();
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() &&
        {
            unreachable();
        }
        
        
        
        int getdimin() const;

        int getdimout() const;
        

        virtual void check() const;

        virtual void print( std::ostream& os ) const;

        bool issquare() const;
        
        /* Apply the operator */
        
        /* x := s A y */
        virtual FloatVector apply( const FloatVector& src, Float scaling = 1. ) const;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling = 1. ) const = 0;
        
        FloatVector createinputvector() const;
        FloatVector createoutputvector() const;
        
        
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
