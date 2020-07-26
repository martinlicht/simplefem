#ifndef INCLUDEGUARD_SOLVER_MinimumResidualMethod_METHOD
#define INCLUDEGUARD_SOLVER_MinimumResidualMethod_METHOD


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

                void iterationStart( 
                        FloatVector& x, const FloatVector& b, 
                        FloatVector& r, 
                        FloatVector& p0, FloatVector& s0, 
                        FloatVector& p1, FloatVector& s1,
                        Float& rr
                        ) const;

                void iterationStep( 
                        FloatVector& x, 
                        FloatVector& r, 
                        FloatVector& p0, FloatVector& s0, 
                        FloatVector& p1, FloatVector& s1,
                        FloatVector& p2, FloatVector& s2,
                        Float& rr
                        ) const;

        private: 

                const LinearOperator& A;   
};
  
  
  
  
#endif
  
  
