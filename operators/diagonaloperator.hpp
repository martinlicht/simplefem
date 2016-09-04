#ifndef INCLUDEGUARD_DIAGONALOPERATOR
#define INCLUDEGUARD_DIAGONALOPERATOR

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
public LinearOperator /* every matrix is a linear operator */
{

    public:

        explicit DiagonalOperator( int, const FloatVector& dia );
        explicit DiagonalOperator( int, const ScalingOperator& scaling );
        virtual ~DiagonalOperator();

        virtual void check() const override;
        virtual void print( std::ostream& ) const override;

        FloatVector& getdiagonal();
        const FloatVector& getdiagonal() const;
        
        virtual void applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const override;

    private:

        FloatVector diagonal;
            
};
  
  

inline DiagonalOperator operator*( const DiagonalOperator& left, const DiagonalOperator& right )
{
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