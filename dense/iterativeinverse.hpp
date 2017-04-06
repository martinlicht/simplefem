#ifndef INCLUDEGUARD_ITERATIVEINVERSE
#define INCLUDEGUARD_ITERATIVEINVERSE

#include <vector>

#include "../basic.hpp"

#include "densematrix.hpp"

void newtoniteration( DenseMatrix& out, const DenseMatrix& src, int n );

#endif