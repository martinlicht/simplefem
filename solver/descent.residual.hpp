#ifndef INCLUDEGUARD_RESIDUAL_DESCENT_METHOD
#define INCLUDEGUARD_RESIDUAL_DESCENT_METHOD


#include <iostream>

#include "../basic.hpp"
#include "iterativesolver.hpp"





/************************
****
****  Class for Conjugate Residual Method
****  - instantiates IterativeSolver
****  - features iteration start and iteration step,
****    which can be called as such, or from solve().
****  
************************/




class ResidualDescentMethod
: public IterativeSolver
{

public:
    
    virtual void check() const override;
    virtual void print( std::ostream& ) const override;
    
    virtual void solve( FloatVector&, const FloatVector& ) const override;
    
    void iterationStart( const FloatVector& x, const FloatVector& b, 
                         FloatVector& r, Float& r_norm ) const;

    void iterationStep( FloatVector& x,
                        FloatVector& r, Float& r_norm ) const;
    
    explicit ResidualDescentMethod( const LinearOperator& op );
    virtual ~ResidualDescentMethod();

    private:
    
            int dimension;

};
  
  
  
  
#endif
  
  
  
  