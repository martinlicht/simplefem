#ifndef INCLUDEGUARD_OPERATOR_SCALING
#define INCLUDEGUARD_OPERATOR_SCALING

#include "../basic.hpp"
#include "linearoperator.hpp"


/************************
****
****  Class for Scalings 
****  - instantiates LinearOperator
****  
************************/


class ScalingOperator:
public LinearOperator /* every scaling operation is a linear operator */
{

    public:

        explicit ScalingOperator( int, Float s );
        virtual ~ScalingOperator();

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<ScalingOperator> clone = std::make_shared<ScalingOperator>( *this );
            return clone;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<ScalingOperator> heir = std::make_unique<ScalingOperator>( *this );
            return heir;
        }
        

        virtual void check() const override;
        virtual void print( std::ostream& ) const override;

        Float getscaling() const;
        void setscaling( Float s );

        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;

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