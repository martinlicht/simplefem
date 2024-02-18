#ifndef INCLUDEGUARD_OPERATOR_BLOCKSTRUCTURE
#define INCLUDEGUARD_OPERATOR_BLOCKSTRUCTURE

#include "../basic.hpp"
#include "linearoperator.hpp"


/************************
****
****  Class for Scalings 
****  - instantiates LinearOperator
****  
************************/


class BlockOperator:
public LinearOperator /* every matrix is a linear operator */
{

    public:
        
        explicit BlockOperator( int dimout, int dimin, const std::vector<std::vector<LinearOperator*>>& ops );
        
        virtual ~BlockOperator();
        
        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;

    private:
        
        std::vector<std::vector<LinearOperator*>> ops;
        
        std::vector<int> dimouts;
        std::vector<int>  dimins;
        
        
};
  


  
  

#endif
