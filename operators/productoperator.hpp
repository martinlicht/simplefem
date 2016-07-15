#ifndef INCLUDEGUARD_PRODUCTOPERATOR
#define INCLUDEGUARD_PRODUCTOPERATOR

#include "../basic.hpp"
#include "linearoperator.hpp"
#include "densematrix.hpp"



class ProductOperator:
public LinearOperator /* every matrix is a linear operator */
{

    public:

        explicit ProductOperator( const LinearOperator&, const LinearOperator& );
        virtual ~ProductOperator();

        virtual void check() const override;
        virtual void print( std::ostream& ) const override;

        virtual void applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const override;

    private:

        const LinearOperator& left;
        const LinearOperator& right;
    
};


inline ProductOperator operator*( const LinearOperator& left, const LinearOperator& right )
{
    return ProductOperator( left, right );
}





  
  

#endif