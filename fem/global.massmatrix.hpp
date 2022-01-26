#ifndef INCLUDEGUARD_FEM_FEECBROKENMASSMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENMASSMATRIX


// #include <cassert>
#include <iostream>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/cholesky.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/matrixtensorproduct.hpp"
#include "../operators/linearoperator.hpp"
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

SparseMatrix FEECBrokenMassMatrix( const Mesh& mesh, int n, int k, int r );

SparseMatrix FEECBrokenMassMatrixRightFactor( const Mesh& mesh, int n, int k, int r );

FloatVector FEECBrokenMassMatrix_cellwisemass( const Mesh& mesh, int n, int k, int r, const FloatVector vec );

SparseMatrix FEECBrokenMassMatrix_cellwiseinverse( const Mesh& mesh, int n, int k, int r );









// inline DenseMatrix elementmassmatrix( int n, int ambientdim, int r, int k, DenseMatrix Jacobian )
// {
//     
//     assert( Jacobian.getdimin() == n && Jacobian.getdimout() == ambientdim ); 
//     assert( n <= ambientdim );
//     assert( n >= 0 && n >= k );
//         
//     DenseMatrix polyMM = polynomialmassmatrix( n, r );
//     
//     DenseMatrix formMM = SubdeterminantMatrix( mesh.getGradientProductMatrix( n, t ), k );
//         
//     return TensorProduct( polyMM, formMM ) * absolute( determinant( Jacobian ) ) / factorial_numerical( n );
//     
// }



#endif
