#ifndef INCLUDEGUARD_DENSE_ITERATIVEINVERSE
#define INCLUDEGUARD_DENSE_ITERATIVEINVERSE

#include <vector>

#include "../basic.hpp"

#include "densematrix.hpp"

void newtoniteration( DenseMatrix& out, const DenseMatrix& src, int n );

#endif