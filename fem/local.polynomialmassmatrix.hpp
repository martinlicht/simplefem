#ifndef INCLUDEGUARD_FEM_POLYNOMIALMASSMATRIX_HPP
#define INCLUDEGUARD_FEM_POLYNOMIALMASSMATRIX_HPP


// #include <cassert>
#include <iostream>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"
#include "../operators/linearoperator.hpp"



DenseMatrix polynomialmassmatrix( int n, int r );


#endif
