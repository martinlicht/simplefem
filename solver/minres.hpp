#ifndef INCLUDEGUARD_SOLVER_MinimumResidualMethod_METHOD
#define INCLUDEGUARD_SOLVER_MinimumResidualMethod_METHOD


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




class MinimumResidualMethod
: public IterativeSolver
{

        public:
        
                explicit MinimumResidualMethod( const LinearOperator& op );
                virtual ~MinimumResidualMethod();

                virtual void check() const override;
                virtual void print( std::ostream& ) const override;
                
                virtual void solve( FloatVector&, const FloatVector& ) const override;

        private: 

                const LinearOperator& A;   
};
  
  
  
  
#endif
  
  
