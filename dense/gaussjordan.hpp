#ifndef INCLUDEGUARD_DENSE_GAUSSJORDAN_HPP
#define INCLUDEGUARD_DENSE_GAUSSJORDAN_HPP

#include <cassert>

#include "../basic.hpp"

#include "densematrix.hpp"


DenseMatrix GaussJordan( DenseMatrix mat );

DenseMatrix GaussJordanInplace( DenseMatrix mat, bool pivoting = true );


 

#endif
