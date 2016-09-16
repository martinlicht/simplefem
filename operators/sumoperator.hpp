#ifndef INCLUDEGUARD_SUMOPERATOR
#define INCLUDEGUARD_SUMOPERATOR

#include "../basic.hpp"
#include "linearoperator.hpp"
#include "densematrix.hpp"


/************************
****
****  Class for Arithmetic Sum of Operators  
****  - instantiates LinearOperator
****  
************************/



class SumOperator:
public LinearOperator /* every matrix is a linear operator */
{

    public:

        explicit SumOperator( const LinearOperator&, const LinearOperator& );
        explicit SumOperator( const LinearOperator&, const LinearOperator&, Float, Float );
        virtual ~SumOperator();

        virtual void check() const override;
        virtual void print( std::ostream& ) const override;

        virtual FloatVector apply( const FloatVector& src, Float scaling ) const override;

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