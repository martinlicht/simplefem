#ifndef INCLUDEGUARD_SOLVER_ITERATIVESOLVER
#define INCLUDEGUARD_SOLVER_ITERATIVESOLVER


#include <iostream>

#include "../basic.hpp"
#include "../operators/linearoperator.hpp"


/************************
****
****  Abstract class for iterative solvers  
****  - uses iteration counter, error tolerance, and internal residual vector 
****  
************************/


  
class IterativeSolver
{

    protected:

        const LinearOperator& internalOperator; 
        mutable FloatVector residual;

    public:  

        explicit IterativeSolver( const LinearOperator& );
        virtual ~IterativeSolver();

        virtual void check() const;
        virtual void print( std::ostream& ) const;

        const LinearOperator& getInternalOperator() const;
        const FloatVector& getResidualVector() const;

        virtual void solve( FloatVector& unknown, const FloatVector& rhs ) const = 0;

        mutable Float tolerance;
        mutable Float recent_deviation;
        mutable int max_iteration_count;
        mutable int recent_iteration_count;
        mutable int print_modulo;

  };
  
 
  
  
  
  
  
#endif
  
  
  
  
