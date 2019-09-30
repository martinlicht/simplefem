#ifndef INCLUDEGUARD_SOLVER_PRECON_CONJUGATERESIDUAL_METHOD
#define INCLUDEGUARD_SOLVER_PRECON_CONJUGATERESIDUAL_METHOD


#include <iostream>

#include "../basic.hpp"
#include "../operators/simpleoperators.hpp"
#include "iterativesolver.hpp"




/************************
****
****  Class for Conjugate Residual Method
****  - instantiates IterativeSolver
****  - features iteration start and iteration step,
****    which can be called as such, or from solve().
****  
************************/




class PreconditionedConjugateResidualMethod
: public IterativeSolver
{

    public:
        
        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        
        virtual void solve( FloatVector&, const FloatVector& ) const override;
        
        void iterationStart( const FloatVector& x, const FloatVector& b, 
                            FloatVector& r, FloatVector& p, FloatVector& Mr, FloatVector& Mp, FloatVector& AMr, FloatVector& AMp,
                            Float& r_MAMnorm ) const;

        void iterationStep( FloatVector& x, const FloatVector& b, 
                            FloatVector& r, FloatVector& p, FloatVector& Mr, FloatVector& Mp, FloatVector& AMr, FloatVector& AMp,
                            Float& r_MAMnorm,
                            FloatVector& MAMp, FloatVector& AMAMp ) const;
        
        explicit PreconditionedConjugateResidualMethod( const LinearOperator& op, const LinearOperator& invprecon );
        virtual ~PreconditionedConjugateResidualMethod();

    private:
    
            int dimension;
            const LinearOperator& invprecon;

};
  
  
  
  
#endif
  
  
  
  
