#ifndef INCLUDEGUARD_MATRIXALGORITHM
#define INCLUDEGUARD_MATRIXALGORITHM

#include <vector>

#include "../basic.hpp"

#include "densematrix.hpp"


/***********************
****  
****  Algorithm for nonlinear dense matrix operations 
****  
****  - Inverse and Determinant 
****  - Polar Decomposition 
****  - Cholesky 
****  - Determinant of Upper triangular 
****  
***********************/


void InverseAndDeterminant( const DenseMatrix&, DenseMatrix&, Float& );

void PolarDecomposition( const DenseMatrix&, DenseMatrix&, DenseMatrix& );

void PolarDecompositionRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );

DenseMatrix CholeskyDecomposition( const DenseMatrix& src );

DenseMatrix UpperTriangularInverse( const DenseMatrix& );        

Float UpperTriangularDeterminant( const DenseMatrix& A );

#endif