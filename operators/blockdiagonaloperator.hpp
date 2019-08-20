#ifndef INCLUDEGUARD_OPERATOR_BLOCKDIAGONAL
#define INCLUDEGUARD_OPERATOR_BLOCKDIAGONAL

#include "../basic.hpp"
#include "linearoperator.hpp"


/************************
****
****  Class for Scalings 
****  - instantiates LinearOperator
****  
************************/


class BlockDiagonalOperator:
public LinearOperator /* every matrix is a linear operator */
{

    public:
        
        explicit BlockDiagonalOperator( int dimout, int dimin, const std::vector<LinearOperator*>& ops );
        virtual ~BlockDiagonalOperator();
        
        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        
        virtual FloatVector apply( const FloatVector& src, Float scaling ) const override;
        
    private:
        
        std::vector<LinearOperator*> ops;
        
};
  


  
  

#endif
