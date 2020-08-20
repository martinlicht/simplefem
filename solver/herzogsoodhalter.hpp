#ifndef INCLUDEGUARD_SOLVER_HERZOGSOODHALTER_METHOD
#define INCLUDEGUARD_SOLVER_HERZOGSOODHALTER_METHOD


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




class HerzogSoodhalterMethod
: public IterativeSolver
{

        public:
        
                explicit HerzogSoodhalterMethod( const LinearOperator& op );
                virtual ~HerzogSoodhalterMethod();

                virtual void check() const override;
                virtual void print( std::ostream& ) const override;
                
                virtual void solve( FloatVector&, const FloatVector& ) const override;

        private: 

                const LinearOperator& A;   
};
  
  
  
  
#endif
  
  
  
  
