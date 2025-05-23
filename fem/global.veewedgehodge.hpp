#ifndef INCLUDEGUARD_FEM_FEECBROKENCONTRACTIONMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENCONTRACTIONMATRIX


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

SparseMatrix FEECBrokenVeeMatrix( const Mesh& mesh, int n, int k, int r, int l, int s, const FloatVector& field );

SparseMatrix FEECBrokenWedgeMatrix( const Mesh& mesh, int n, int k, int r, int l, int s, const FloatVector& field );

SparseMatrix FEECBrokenHodgeStarMatrix( const Mesh& mesh, int n, int k, int r );

FloatVector FEECVolumeFormIntegral( const Mesh& mesh, int n, int r );

FloatVector FEECScalarIntegral( const Mesh& mesh, int n, int r );

DenseMatrix EuclideanHodgeStar( int n, int k );

SparseMatrix FEECBrokenHodgeStarMatrix_Alternative( const Mesh& mesh, int n, int k, int r );



#endif
