#ifndef INCLUDEGUARD_SIMPLESOLVER
#define INCLUDEGUARD_SIMPLESOLVER

#include <vector>

#include "../basic.hpp"

#include "densematrix.hpp"

 

DenseMatrix DiagonalPart( const DenseMatrix& A );
DenseMatrix DiagonalInverse( const DenseMatrix& A );
DenseMatrix InvertDiagonal( DenseMatrix& A );
Float DiagonalDeterminant( const DenseMatrix& A );
void DiagonalSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );

DenseMatrix LowerTriangularPart( const DenseMatrix& A );
// DenseMatrix LowerTriangularInverse( const DenseMatrix& A );
// DenseMatrix InvertLowerTriangular( const DenseMatrix& A );
Float LowerTriangularDeterminant( const DenseMatrix& A );
// void LowerTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );

DenseMatrix UpperTriangularPart( const DenseMatrix& A );
DenseMatrix UpperTriangularInverse( const DenseMatrix& A );
// DenseMatrix InvertUpperTriangular( const DenseMatrix& A );
Float UpperTriangularDeterminant( const DenseMatrix& A );
// void UpperTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );

DenseMatrix LowerUnitTriangularPart( const DenseMatrix& A );
// DenseMatrix LowerUnitTriangularInverse( const DenseMatrix& A );
// void LowerUnitTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );

DenseMatrix UpperUnitTriangularPart( const DenseMatrix& A );
// DenseMatrix UpperUnitTriangularInverse( const DenseMatrix& A );
// void UpperUnitTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );



#endif