#ifndef INCLUDEGUARD_SOLVER_CONJUGATEGRADIENT_METHOD
#define INCLUDEGUARD_SOLVER_CONJUGATEGRADIENT_METHOD


#include <ostream>

#include "../basic.hpp"
#include "../operators/linearoperator.hpp"
#include "iterativesolver.hpp"




/************************
****
****  Class for Conjugate Gradient Method
****  - instantiates IterativeSolver
****  - features iteration start and iteration step,
****    which can be called as such, or from solve().
****  
************************/




class ConjugateGradientMethod
: public IterativeSolver
{

        public:
        
                explicit ConjugateGradientMethod( const LinearOperator& op );
                virtual ~ConjugateGradientMethod();

                virtual void check() const override;
                virtual void print( std::ostream& ) const override;
                
                virtual void solve( FloatVector&, const FloatVector& ) const override;

        private: 

                const LinearOperator& A;   
};
  
  
  
  
#endif
  
  
  
  
