#ifndef INCLUDEGUARD_QR_FACTORIZATION
#define INCLUDEGUARD_QR_FACTORIZATION

#include <vector>

#include "../basic.hpp"

#include "densematrix.hpp"


 // QR decomposition
 // LQ decomposition 
 // QR repeated 
 // LQ repeated 
 
void QRFactorization( const DenseMatrix&, DenseMatrix&, DenseMatrix& );

void QRFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );

void LQFactorization( const DenseMatrix&, DenseMatrix&, DenseMatrix& );

void LQFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );



#endif