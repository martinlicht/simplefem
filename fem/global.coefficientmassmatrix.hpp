#ifndef INCLUDEGUARD_FEM_FEECBROKENMASSMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENMASSMATRIX


#include <functional>

#include "../operators/floatvector.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/local.polynomialmassmatrix.hpp"



//////////////////////////////////////////////////////
//                                                  //
//  Matrix for the mass pairing                     //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECBrokenCoefficientMassMatrix( const Mesh& mesh, int n, int k, int r, int s, std::function<DenseMatrix(FloatVector)>& generator );

#endif
