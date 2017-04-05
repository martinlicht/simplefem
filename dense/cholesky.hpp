#ifndef INCLUDEGUARD_DENSE_CHOLESKY
#define INCLUDEGUARD_DENSE_CHOLESKY

#include <vector>

#include "../basic.hpp"

#include "densematrix.hpp"


DenseMatrix CholeskyDecomposition( const DenseMatrix& src );

DenseMatrix CholeskyDecompositionBanachchiewicz( const DenseMatrix& src );

// TODO: Cholesky with Crout pattern, and other possible patterns 
// TODO: Break down condition?

// TODO: Cholesky with Pivoting 


#endif