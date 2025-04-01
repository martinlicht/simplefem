#ifndef INCLUDEGUARD_FEM_CONFORMING
#define INCLUDEGUARD_FEM_CONFORMING

#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"


//////////////////////////////////////////////////////
//                                                  //
//  Matrix for inclusion of                         //
//  continuous Sullivan forms                       //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECSullivanInclusionMatrix( const Mesh& mesh, int n, int k, int r );

//////////////////////////////////////////////////////
//                                                  //
//  Matrix for inclusion of                         //
//  continuous Whitney forms                        //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECWhitneyInclusionMatrix( const Mesh& mesh, int n, int k, int r );



//////////////////////////////////////////////////////
//                                                  //
//  Matrix for inclusion of                         //
//  continuous AFW forms                            //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECAFWInclusionMatrix( const Mesh& mesh, int n, int k, int r );




//////////////////////////////////////////////////////
//                                                  //
//  Matrix for inclusion of                         //
//  continuous Sullivan forms                       //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

DiagonalOperator FEECSullivanFlagMatrix( const Mesh& mesh, int n, int k, int r );

DiagonalOperator FEECWhitneyFlagMatrix( const Mesh& mesh, int n, int k, int r );








//////////////////////////////////////////////////////
//                                                  //
//  Matrix for inclusion of                         //
//  continuous Sullivan forms                       //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

enum class FEECAveragingMode : uint8_t
{
    weighted_uniformly,
    arbitrary_choice,
    weighted_by_volume
};

SparseMatrix FEECSullivanAveragingMatrix( const Mesh& mesh, int n, int k, int r, FEECAveragingMode mode );








#endif
