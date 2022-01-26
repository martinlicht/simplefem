#ifndef INCLUDEGUARD_FEM_FEECBROKENDIFFMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENDIFFMATRIX


// #include <cassert>
#include <iostream>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/matrixtensorproduct.hpp"
#include "../operators/linearoperator.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/local.polynomialmassmatrix.hpp"


//////////////////////////////////////////////////////
//                                                  //
//  Matrix for the exterior derivative              //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECBrokenDiffMatrix( const Mesh& mesh, int n, int k, int r );




#endif
