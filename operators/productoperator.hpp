#ifndef INCLUDEGUARD_OPERATOR_PRODUCT
#define INCLUDEGUARD_OPERATOR_PRODUCT

#include "../basic.hpp"
#include "linearoperator.hpp"
#include "scalingoperator.hpp"


/************************
****
****  Class for Products of LinearOperators 
****  - instantiates LinearOperator
****  
************************/

class ProductOperator:
public LinearOperator /* every matrix is a linear operator */
{

    public:

        explicit ProductOperator( const LinearOperator&, const LinearOperator& );
        virtual ~ProductOperator();

        virtual void check() const override;
        virtual void print( std::ostream& ) const override;

        virtual FloatVector apply( const FloatVector& src, Float scaling = 1. ) const override;

    private:

        const LinearOperator& left;
        const LinearOperator& right;
    
};


inline ProductOperator operator*( const LinearOperator& left, const LinearOperator& right )
{
    left.check(); right.check();
    return ProductOperator( left, right );
}

inline ProductOperator operator*( Float left, const LinearOperator& right )
{
    return ProductOperator( ScalingOperator(right.getdimout(),left), right );
}

inline ProductOperator operator*( const LinearOperator& left, Float right )
{
    return ProductOperator( left, ScalingOperator(left.getdimin(),right) );
}



  
  

#endif
