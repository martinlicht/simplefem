#ifndef INCLUDEGUARD_FEM_LAGRANGE_MATRICES
#define INCLUDEGUARD_FEM_LAGRANGE_MATRICES

#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"






//////////////////////////////////////////////////////
//                                                  //
//  Matrix for the broken Lagrange mass pairing     //
//  with poly degree r                              //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix LagrangeBrokenMassMatrix( const Mesh& mesh, int r );



///////////////////////////////////////////////////////
//                                                   //
//  Matrix for the continuous Lagrange mass pairing  //
//  with poly degree r                               //
//                                                   //
///////////////////////////////////////////////////////

SparseMatrix LagrangeMassMatrix( const Mesh& mesh, int r );



//////////////////////////////////////////////////////////
//                                                      //
//  Matrix for the broken Lagrange stiffness pairing    //
//  with poly degree r                                  //
//                                                      //
//////////////////////////////////////////////////////////

SparseMatrix LagrangeBrokenStiffnessMatrix( const Mesh& mesh, int r );






//////////////////////////////////////////////////////////////
//                                                          //
//  Matrix for the continuous Lagrange stiffness pairing    //
//  with poly degree r                                      //
//                                                          //
//////////////////////////////////////////////////////////////

SparseMatrix LagrangeStiffnessMatrix( const Mesh& mesh, int r );






////////////////////////////////////////////////
//                                            //
//  Matrix for the inclusion of               //
//  continuous Lagrange elements              //
//  into the broken space                     //
//  with poly degree r                        //
//                                            //
////////////////////////////////////////////////

SparseMatrix LagrangeInclusionMatrix( const Mesh& mesh, int n, int r );







#endif
