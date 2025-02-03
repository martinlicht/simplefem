#ifndef INCLUDEGUARD_FEM_FEECBROKENHODGESTARMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENHODGESTARMATRIX

#include "../dense/densematrix.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"


//////////////////////////////////////////////////////
//                                                  //
//  Matrix for the mass pairing                     //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

DenseMatrix EuclideanHodgeStar( int n, int k );

SparseMatrix FEECBrokenHodgeStarPairingMatrix( const Mesh& mesh, int n, int k, int r );

#endif
