#ifndef INCLUDEGUARD_FEM_UNPHYSICAL_HPP
#define INCLUDEGUARD_FEM_UNPHYSICAL_HPP

#include <algorithm>
#include <vector>

#include "../basic.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"


/************************
****
****  Unphysical operations
****  
************************/


SparseMatrix FEECCanonicalizeBroken( const Mesh& mesh, int n, int k, int r );

SparseMatrix FEECRandomizeBroken( const Mesh& mesh, int n, int k, int r );
 

#endif
