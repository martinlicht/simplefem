#ifndef INCLUDEGUARD_DENSE_FUNCTIONS_HPP
#define INCLUDEGUARD_DENSE_FUNCTIONS_HPP

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

DenseMatrix Inverse( DenseMatrix );

void Inverse_InSitu( DenseMatrix& );

// void InverseAndDeterminant( const DenseMatrix&, DenseMatrix&, Float& );

DenseMatrix SubdeterminantMatrixSquare( const DenseMatrix& A, int k );

DenseMatrix SubdeterminantMatrix( const DenseMatrix& A, int k );



/*
 *
 * Internal functions 
 * 
 */

inline Float Determinant_laplaceexpansion( const DenseMatrix& );

inline Float Determinant_gauss( DenseMatrix );

inline void Inverse_CramersRule_InSitu( DenseMatrix& );

inline void Inverse_gauss_InSitu( DenseMatrix&, bool pivoting = true );










#endif
