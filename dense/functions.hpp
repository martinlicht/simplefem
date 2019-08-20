#ifndef INCLUDEGUARD_DENSE_FUNCTIONS
#define INCLUDEGUARD_DENSE_FUNCTIONS

#include <vector>

#include "../basic.hpp"

#include "densematrix.hpp"

 
 // Tranpose, giving new matrix 
 // Tranpose, in-place

DenseMatrix Transpose( const DenseMatrix& ); 
void TransposeInSitu( DenseMatrix& );
DenseMatrix TransposeSquare( const DenseMatrix& ); 
void TransposeSquareInSitu( DenseMatrix& ); 








 // determinant 
 // cofactor matrix 
 // Inverse 
 // Inverse and determinant 
 // SubdeterminantMatrix

Float Determinant( const DenseMatrix& );

DenseMatrix CofactorMatrix( const DenseMatrix& ); 

DenseMatrix Inverse( const DenseMatrix& );

void InverseAndDeterminant( const DenseMatrix&, DenseMatrix&, Float& );

DenseMatrix SubdeterminantMatrix( const DenseMatrix& A, int k );



/*
 *
 * Internal functions 
 * 
 */

Float Determinant_laplaceexpansion( const DenseMatrix& );










#endif