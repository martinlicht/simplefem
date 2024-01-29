#ifndef INCLUDEGUARD_FEM_TRACEMATRIX
#define INCLUDEGUARD_FEM_TRACEMATRIX


#include <vector>
#include <iostream>
// #include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"


SparseMatrix BrokenTraceMatrix( const Mesh& mesh, int n, int k, int r, bool is_signed );



#endif
