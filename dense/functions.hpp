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

// // // // // // // // matrix trace 
// // // // // // // 
// // // // // // // Float MatrixTrace( const DenseMatrix& x );
// // // // // // // 
// // // // // // // // Gerschgorin circles : row/column 
// // // // // // // 
// // // // // // // DenseMatrix Gerschgorin( const DenseMatrix& );
// // // // // // // DenseMatrix GerschgorinRow( const DenseMatrix& );
// // // // // // // DenseMatrix GerschgorinColumn( const DenseMatrix& );
// // // // // // // 
// // // // // // // // Matrix norms 
// // // // // // // 
// // // // // // // Float NormL1( const DenseMatrix& );
// // // // // // // Float NormFrobenius( const DenseMatrix& );
// // // // // // // Float NormMax( const DenseMatrix& );
// // // // // // // Float NormLp( const DenseMatrix&, Float p );
// // // // // // // Float NormRowCol( const DenseMatrix&, Float p, Float q );
// // // // // // // Float NormColRow( const DenseMatrix&, Float p, Float q );
// // // // // // // 
// // // // // // // Float NormOperatorL1( const DenseMatrix& );
// // // // // // // Float NormOperatorMax( const DenseMatrix& );







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