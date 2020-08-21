#ifndef INCLUDEGUARD_SOLVER_PRECON_CONJUGATERESIDUAL_METHOD
#define INCLUDEGUARD_SOLVER_PRECON_CONJUGATERESIDUAL_METHOD


#include <ostream>

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
        
        explicit PreconditionedConjugateResidualMethod( const LinearOperator& op, const LinearOperator& M );
        virtual ~PreconditionedConjugateResidualMethod();

        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        
        virtual void solve( FloatVector&, const FloatVector& ) const override;
        
    private:

        const LinearOperator& A;   
        const LinearOperator& M;

};
  
  
  
  
#endif
  
  
  
  
