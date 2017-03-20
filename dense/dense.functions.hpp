#ifndef INCLUDEGUARD_DENSE_FUNCTIONS
#define INCLUDEGUARD_DENSE_FUNCTIONS

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
 
 // Tranpose, giving new matrix 
 // Tranpose, in-place

 DenseMatrix Transpose( const DenseMatrix& ); 
 void TransposeInSitu( DenseMatrix& );
 DenseMatrix TransposeSquare( const DenseMatrix& ); 
 void TransposeSquareinSitu( DenseMatrix& ); 

 
 // matrix trace 
 
 // Gerschgorin circles : row/column 
 
 // Norms: Frobenius, row norm, column norm,
 
 // Diagonal inverse
 // Left triangular inverse 
 // Right triangular inverse 
 // Unit Left triangular inverse
 // Unit Right triangular inverse 
 
 // Diagonal solve 
 // Left triangular solve 
 // Right triangular solve 
 // Unit Left triangular solve 
 // Unit Right triangular solve 
 
 // (Weighting between left and right triangular solves)
  
 
 // determinant 
 // adjunct matrix 
 // Inverse 
 // Solve 
 // Inverse and determinant 
 // Tensorproduct matrix 
 // SubdeterminantMatrix



        
        


void InverseAndDeterminant( const DenseMatrix&, DenseMatrix&, Float& );

DenseMatrix Inverse( const DenseMatrix& );

Float Determinant( const DenseMatrix& );

DenseMatrix UpperTriangularInverse( const DenseMatrix& A );        

Float UpperTriangularDeterminant( const DenseMatrix& A );

DenseMatrix MatrixTensorProduct( const DenseMatrix& left, const DenseMatrix& right );

DenseMatrix SubdeterminantMatrix( const DenseMatrix& A, int k );

DenseMatrix adjunctMatrix( const DenseMatrix& ); 




#endif