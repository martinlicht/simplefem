#ifndef INCLUDEGUARD_OPERATOR_SUM
#define INCLUDEGUARD_OPERATOR_SUM

#include "../basic.hpp"
#include "linearoperator.hpp"


/************************
****
****  Class for Arithmetic Sum of Operators  
****  - instantiates LinearOperator
****  
************************/



class SumOperator:
public LinearOperator 
{

    public:

        explicit SumOperator( const LinearOperator& left, const LinearOperator& right, Float leftscale = 1., Float rightscale = 1. );
        virtual ~SumOperator();

        virtual void check() const override;
        virtual void print( std::ostream& ) const override;

        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;

    private:

        const LinearOperator& left;
        const LinearOperator& right;
        Float leftscale, rightscale;
    
};


inline SumOperator operator+( const LinearOperator& left, const LinearOperator& right )
{
    left.check(); right.check();
    return SumOperator( left, right );
}





  
  

#endif
