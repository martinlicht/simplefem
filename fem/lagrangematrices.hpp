#ifndef INCLUDEGUARD_FEM_LAGRANGE_MATRICES
#define INCLUDEGUARD_FEM_LAGRANGE_MATRICES



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
#include "../operators/linearoperator.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../dense/matrixtensorproduct.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/local.polynomialmassmatrix.hpp"







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
