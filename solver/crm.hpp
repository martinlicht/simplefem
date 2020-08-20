#ifndef INCLUDEGUARD_SOLVER_CONJUGATERESIDUAL_METHOD
#define INCLUDEGUARD_SOLVER_CONJUGATERESIDUAL_METHOD


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




class ConjugateResidualMethod
: public IterativeSolver
{

        public:
        
                explicit ConjugateResidualMethod( const LinearOperator& op );
                virtual ~ConjugateResidualMethod();

                virtual void check() const override;
                virtual void print( std::ostream& ) const override;
                
                virtual void solve( FloatVector&, const FloatVector& ) const override;

                virtual void solve_robust( FloatVector&, const FloatVector& ) const;

                virtual void solve_fast( FloatVector&, const FloatVector& ) const;

        private: 

                const LinearOperator& A;   
};
  
  
  
  
#endif
  
  
  
  
