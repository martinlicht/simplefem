#ifndef INCLUDEGUARD_FEM_FEECMETRIC
#define INCLUDEGUARD_FEM_FEECMETRIC


#include "../operators/floatvector.hpp"
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

SparseMatrix FEECBrokenMassMatrix( const Mesh& mesh, int n, int k, int r );

SparseMatrix FEECBrokenMassMatrixRightFactor( const Mesh& mesh, int n, int k, int r );

FloatVector FEECBrokenMassMatrix_cellwisemass( const Mesh& mesh, int n, int k, int r, const FloatVector& vec );

SparseMatrix FEECBrokenMassMatrix_cellwiseinverse( const Mesh& mesh, int n, int k, int r );


















//////////////////////////////////////////////////////
//                                                  //
//  Matrix for the mass pairing                     //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECBrokenCoefficientMassMatrix( const Mesh& mesh, int n, int k, int r,
                                              int w, const std::function<DenseMatrix(const FloatVector&)>& generator );













//////////////////////////////////////////////////////
//                                                  //
//  Matrix for the mass pairing                     //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECBrokenVeeMatrix( const Mesh& mesh, int n, int k, int r, int l, int s, const FloatVector& field );

SparseMatrix FEECBrokenWedgeMatrix( const Mesh& mesh, int n, int k, int r, int l, int s, const FloatVector& field );

SparseMatrix FEECBrokenHodgeStarMatrix( const Mesh& mesh, int n, int k, int r );

FloatVector FEECVolumeFormIntegral( const Mesh& mesh, int n, int r );

FloatVector FEECScalarIntegral( const Mesh& mesh, int n, int r );

DenseMatrix EuclideanHodgeStar( int n, int k );

SparseMatrix FEECBrokenHodgeStarMatrix_Alternative( const Mesh& mesh, int n, int k, int r );



#endif
