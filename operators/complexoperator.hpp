#ifndef INCLUDEGUARD_OPERATOR_COMPLEXOPERATOR_HPP
#define INCLUDEGUARD_OPERATOR_COMPLEXOPERATOR_HPP

#include <cassert>
#include <functional>
#include <memory>
#include <ostream>
#include <utility>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


/************************
****
****  Class for Complex Operators 
****  
************************/


FloatVector RealPart( const FloatVector& vec );
FloatVector ImaginaryPart( const FloatVector& vec );
FloatVector ComplexFloatVector( const FloatVector& real, const FloatVector& imag );

class ComplexOperator final
: public LinearOperator 
{

    public:

        /* Constructors */
        
        explicit ComplexOperator( const LinearOperator& real );
        explicit ComplexOperator( const LinearOperator& real, const LinearOperator& imag );

        /* standard methods for operators */
        
        ComplexOperator()                                         = delete;
        ComplexOperator( const ComplexOperator& )                = default;
        ComplexOperator( ComplexOperator&& )                     = default;
        ComplexOperator& operator=( const ComplexOperator& vec ) = default;
        ComplexOperator& operator=( ComplexOperator&& vec )      = default; 

        virtual ~ComplexOperator();
        
        /* standard interface */
        
        virtual void check() const override;
        
        virtual std::string text() const override;
        
        virtual std::string text( const bool embellish ) const;

        /* OTHER METHODS */
        
        virtual ComplexOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        

        const LinearOperator& real() const;
        const LinearOperator& imaginary() const;
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;
        
    private:

        const LinearOperator& part_real;
        const LinearOperator& part_imag;
            
};



  
  
  
  

#endif
