#ifndef INCLUDEGUARD_DENSE_FACTORIZATION
#define INCLUDEGUARD_DENSE_FACTORIZATION

#include <vector>

#include "../basic.hpp"

#include "densematrix.hpp"


/***********************
****  
****  Algorithm for dense matrix operations 
****  
****  - Inverse and Determinant 
****  - QR Decomposition 
****  - Cholesky 
****  - Determinant of Upper triangular 
****  - tensor product matrix 
****  
***********************/
 
 // Stabilized Gram-Schmidt 
 
 // QR decomposition
 // LQ decomposition 
 // QR repeated 
 // LQ repeated 
 
 // LR factorization 
 // LR factorization, row pivot  
 // LR factorization, column pivot  
 // LR factorization, full pivot
 // Dolittle vs Crout  
 
 // Cholesky
 // Cholesky, pivot
 
void QRFactorization( const DenseMatrix&, DenseMatrix&, DenseMatrix& );

void QRFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );

DenseMatrix CholeskyDecomposition( const DenseMatrix& src );



#endif