#ifndef INCLUDEGUARD_MATRIXALGORITHM
#define INCLUDEGUARD_MATRIXALGORITHM

#include <vector>
#include "../basic.hpp"
#include "densematrix.hpp"


// DenseMatrix Inverse( const DenseMatrix& );

void PolarDecomposition( const DenseMatrix&, DenseMatrix&, DenseMatrix& );

void PolarDecompositionRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );

DenseMatrix CholeskyDecomposition( const DenseMatrix& src );

// Float GaussFactorization( const DenseMatrix&, DenseMatrix&, DenseMatrix& );

// Float CholeskyFactorization( const DenseMatrix&, DenseMatrix& );

// DenseMatrix LowerTriangularInverse( const DenseMatrix& );        

DenseMatrix UpperTriangularInverse( const DenseMatrix& );        


#endif