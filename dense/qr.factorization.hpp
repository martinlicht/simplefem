#ifndef INCLUDEGUARD_QR_FACTORIZATION
#define INCLUDEGUARD_QR_FACTORIZATION

#include <vector>

#include "../basic.hpp"

#include "densematrix.hpp"


 // QR decomposition
 // LQ decomposition 
 
void QRFactorization( const DenseMatrix&, DenseMatrix&, DenseMatrix& );

void LQFactorization( const DenseMatrix&, DenseMatrix&, DenseMatrix& );



 // QR repeated 
 // LQ repeated 

// void QRFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );
// void LQFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );



#endif