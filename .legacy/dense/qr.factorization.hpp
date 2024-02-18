#ifndef INCLUDEGUARD_DENSE_QRFACTORIZATION_HPP
#define INCLUDEGUARD_DENSE_QRFACTORIZATION_HPP

#include "../basic.hpp"

#include "densematrix.hpp"


 // QR decomposition
 // LQ decomposition 
 
void QRFactorization( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R );

void LQFactorization( const DenseMatrix& A, DenseMatrix& L, DenseMatrix& Q );

FloatVector SolveOverconstrained( const DenseMatrix& A, const FloatVector& v );


 // QR repeated 
 // LQ repeated 

// void QRFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );
// void LQFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );



#endif
