#ifndef INCLUDEGUARD_FEM_FEECSULLIVANINCLUSIONMATRIX
#define INCLUDEGUARD_FEM_FEECSULLIVANINCLUSIONMATRIX


#include <iostream>
#include <utility>
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

#include "indexfunctions.hpp"



//////////////////////////////////////////////////////
//                                                  //
//  Matrix for inclusion of                         //
//  continuous Sullivan forms                       //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//  and lifts up rplus degrees                      //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECSullivanInclusionMatrix( const Mesh& mesh, int n, int k, int r );




#endif
