#ifndef INCLUDEGUARD_SOLVER_CONJUGATERESIDUAL_METHOD
#define INCLUDEGUARD_SOLVER_CONJUGATERESIDUAL_METHOD


#include <iostream>

#include "../basic.hpp"
#include "../operators/diagonaloperator.hpp"
#include "iterativesolver.hpp"




/************************
****
****  Class for Conjugate Residual Method
****  - instantiates IterativeSolver
****  - features iteration start and iteration step,
****    which can be called as such, or from solve().
****  
************************/




class ConjugateResidualMethod
: public IterativeSolver
{

        public:
        
        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        
        virtual void solve( FloatVector&, const FloatVector& ) const override;
        
        void iterationStart( const FloatVector& x, const FloatVector& b, 
                                FloatVector& r, FloatVector& d, FloatVector& Ar, FloatVector& Ad,
                                Float& rAnorm ) const;

        void iterationStep( FloatVector& x,
                                FloatVector& r, FloatVector& d, FloatVector& Ar, FloatVector& Ad,
                                Float& rAnorm, FloatVector& p ) const;
        
        explicit ConjugateResidualMethod( const LinearOperator& op );
        virtual ~ConjugateResidualMethod();

    private:
    
            int dimension;
    
};
  
  
  
  
#endif
  
  
  
  