#ifndef INCLUDEGUARD_NORMALEQUATIONS_CONJUGATERESIDUAL_METHOD
#define INCLUDEGUARD_NORMALEQUATIONS_CONJUGATERESIDUAL_METHOD


#include <iostream>

#include "../basic.hpp"
#include "iterativesolver.hpp"
#include "densematrix.hpp"
#include "crm.hpp"




/************************
****
****  Class for Conjugate Residual Method on Normal Equations
****  - instantiates IterativeSolver
****  - specialized to Dense Matrices
****  
************************/


// In order to make this class an instance of IterativeSolver,
// we have to deal with the iterativesolver parent class being constructed
// before the member variables. 
// This could be done by several means.
// (a) initiliaze the parent class with a ProductOperator, then retrieve data later
// (b) equip iterativesolver with reinitialize()
// (c) busy evalualtion / bad code 
// (d) delegate constructors 


class NormalEquationsConjugateResidualMethod
{
    
    private:
        
        NormalEquationsConjugateResidualMethod( const DenseMatrix& op, const DenseMatrix& opt )
        : NormalEquationsConjugateResidualMethod( op, opt, opt * op ){ };
        
        NormalEquationsConjugateResidualMethod( const DenseMatrix& op, const DenseMatrix& opt, const DenseMatrix& opsys );
        /* this is the actual constructor */
        
    public:
    
        explicit NormalEquationsConjugateResidualMethod( const DenseMatrix& op )
        : NormalEquationsConjugateResidualMethod( op, op.transpose() ){ };
        
        ~NormalEquationsConjugateResidualMethod();
        
        void check() const;
        void print( std::ostream& ) const;
        
        void solve( FloatVector&, const FloatVector& ) const;
    
    private:
    
        const DenseMatrix& matrix_original;
        const DenseMatrix matrix_transposed;
        const DenseMatrix matrix_system;
        ConjugateResidualMethod crm;

};
  
  
  
  
#endif
  
  
  
  