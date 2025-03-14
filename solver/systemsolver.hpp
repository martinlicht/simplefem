#ifndef INCLUDEGUARD_SOLVER_SYSTEMSOLVER
#define INCLUDEGUARD_SOLVER_SYSTEMSOLVER


#include "../base/include.hpp"
#include "../operators/linearoperator.hpp"

int BlockHerzogSoodhalterMethod( 
    FloatVector& x_A, 
    FloatVector& x_C, 
    const FloatVector& b_A, 
    const FloatVector& b_C, 
    const LinearOperator& A, const LinearOperator& Bt, const LinearOperator& B, const LinearOperator& C, 
    Float precision,
    int print_modulo,
    const LinearOperator& PAinv, const LinearOperator& PCinv
);
 
  
  
  
  
  
#endif
