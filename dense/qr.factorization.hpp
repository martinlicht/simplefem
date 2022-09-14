#ifndef INCLUDEGUARD_DENSE_QRFACTORIZATION_HPP
#define INCLUDEGUARD_DENSE_QRFACTORIZATION_HPP

#include "../basic.hpp"

#include "densematrix.hpp"


 // QR decomposition
 // LQ decomposition 
 
void QRFactorization( const DenseMatrix&, DenseMatrix&, DenseMatrix& );

void LQFactorization( const DenseMatrix&, DenseMatrix&, DenseMatrix& );

FloatVector SolveOverconstrained( const DenseMatrix& A, const FloatVector& v );


 // QR repeated 
 // LQ repeated 

// void QRFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );
// void LQFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );



#endif
