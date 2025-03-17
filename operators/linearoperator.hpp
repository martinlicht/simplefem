#ifndef INCLUDEGUARD_OPERATOR_LINEAROPERATOR_HPP
#define INCLUDEGUARD_OPERATOR_LINEAROPERATOR_HPP




//#include <memory>
// #include <ostream>

#include "../base/include.hpp"
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
        
        /* Constructors */
        
        explicit LinearOperator( int );
        explicit LinearOperator( int, int );

        /* standard methods for operators */
        
        LinearOperator()                             = delete;
        LinearOperator( const LinearOperator& )      = default;
        LinearOperator( LinearOperator&& )           = default;
        LinearOperator& operator=( const LinearOperator& op ) = default;
        LinearOperator& operator=( LinearOperator&& op )      noexcept = default;
        virtual ~LinearOperator()                             noexcept = default;

        /* standard interface */
        
        virtual void check() const;

        virtual std::string text() const = 0;
        
        // void print( std::ostream& os ) const;

        // // void lg() const { LOG << text() << nl; };
        
        /* OTHER METHODS */
        
        virtual LinearOperator* pointer_to_heir() && = 0;        
        
        
        int getdimin() const;

        int getdimout() const;
        
        bool is_square() const;
        
        /* Apply the operator */
        
        /* x := s A y */
        virtual FloatVector apply( const FloatVector& src ) const;
        virtual FloatVector apply( const FloatVector& src, Float scaling ) const;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const = 0;
        
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
    op.apply( ret, vec, 1. );
    vec.check();
    return ret;
}
  
template<typename Stream>
inline decltype(auto) operator<<( Stream&& os, const LinearOperator& op )
{
    op.check();
    os << op.text(); // op.print( os );
    op.check();
    return std::forward<Stream>(os);
}
  


  
  
  
  
#endif
