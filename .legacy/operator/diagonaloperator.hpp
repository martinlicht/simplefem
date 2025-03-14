#ifndef INCLUDEGUARD_OPERATOR_DIAGONAL
#define INCLUDEGUARD_OPERATOR_DIAGONAL

#include "../basic.hpp"
#include "linearoperator.hpp"
#include "scalingoperator.hpp"



/************************
****
****  Class for diagonal matrices 
****  - instantiates LinearOperator
****  - only square matrices 
****  
************************/

class DiagonalOperator:
public LinearOperator 
{

    public:

        explicit DiagonalOperator( int, Float s );
        explicit DiagonalOperator( int, const FloatVector& dia );
        explicit DiagonalOperator( int, const ScalingOperator& scaling );
        virtual ~DiagonalOperator();

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<DiagonalOperator> clone = std::make_shared<DiagonalOperator>( *this );
            return clone;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<DiagonalOperator> heir = std::make_unique<DiagonalOperator>( *this );
            return heir;
        }
        

        virtual void check() const override;
        virtual void print( std::ostream& ) const override;

        FloatVector& getdiagonal();
        const FloatVector& getdiagonal() const;
        
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;

    private:

        FloatVector diagonal;
            
};
  
  

inline DiagonalOperator operator*( const DiagonalOperator& left, const DiagonalOperator& right )
{
    left.check();
    right.check();
    
    assert( left.getdimin() == right.getdimout() );
    assert( left.getdimout() == right.getdimin() );
    
    const FloatVector& leftdia = left.getdiagonal();
    const FloatVector& rightdia = right.getdiagonal();
    const int dimension = leftdia.getdimension();
    
    assert( leftdia.getdimension() == rightdia.getdimension() );
    
    return DiagonalOperator( left.getdimout(), 
                             FloatVector( leftdia.getdimension(), 
                                          [&](int d) -> Float { 
                                            assert( 0 <= d && d < dimension ); 
                                            return leftdia[d] * rightdia[d];
                                          }
                                        )
                           );
}  
  
  
  
  
  
#endif