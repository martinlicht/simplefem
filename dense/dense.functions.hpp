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
void TransposeSquareinSitu( DenseMatrix& ); 

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
 // adjunct matrix 
 // Inverse 
 // Inverse and determinant 
 // SubdeterminantMatrix

void InverseAndDeterminant( const DenseMatrix&, DenseMatrix&, Float& );

DenseMatrix Inverse( const DenseMatrix& );

Float Determinant( const DenseMatrix& );

DenseMatrix adjunctMatrix( const DenseMatrix& ); 

DenseMatrix SubdeterminantMatrix( const DenseMatrix& A, int k );


// Tensorproduct matrix 

DenseMatrix MatrixTensorProduct( const DenseMatrix& left, const DenseMatrix& right );




// TODO: Put triangular solves somewhere else 
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

DenseMatrix UpperTriangularInverse( const DenseMatrix& A );        

Float UpperTriangularDeterminant( const DenseMatrix& A );







#endif