#ifndef INCLUDEGUARD_SOLVER_RESIDUAL_DESCENT_METHOD
#define INCLUDEGUARD_SOLVER_RESIDUAL_DESCENT_METHOD


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
        
                explicit ResidualDescentMethod( const LinearOperator& op );
                virtual ~ResidualDescentMethod();

                virtual void check() const override;
                virtual void print( std::ostream& ) const override;
                
                virtual void solve( FloatVector&, const FloatVector& ) const override;

        private: 

                const LinearOperator& A;   
};


  
  
  
#endif
  
  
  
  