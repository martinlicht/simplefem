#ifndef INCLUDEGUARD_SCALINGOPERATOR
#define INCLUDEGUARD_SCALINGOPERATOR

#include "../basic.hpp"
#include "linearoperator.hpp"


/************************
****
****  Class for Scalings 
****  - instantiates LinearOperator
****  
************************/


class ScalingOperator:
public LinearOperator /* every matrix is a linear operator */
{

    public:

        explicit ScalingOperator( int, Float s );
        virtual ~ScalingOperator();

        virtual void check() const override;
        virtual void print( std::ostream& ) const override;

        Float getscaling() const;
        void setscaling( Float s );

        virtual void applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const override;

    private:

        Float scaling;
    
};
  
  
inline ScalingOperator operator*( const ScalingOperator& left, const ScalingOperator& right )
{
    left.check();
    right.check();
    assert( left.getdimin() == right.getdimout() );
    assert( left.getdimout() == right.getdimin() );
    
    return ScalingOperator( left.getdimout(), left.getscaling() * right.getscaling() );
}  
  
  
  

#endif