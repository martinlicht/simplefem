#ifndef INCLUDEGUARD_FEM_INDEXFUNCTIONS_HPP
#define INCLUDEGUARD_FEM_INDEXFUNCTIONS_HPP


#include <iostream>
#include <utility>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"




//////////////////////////////////////////////////////
//                                                  //
//  List of Sullivan Indices                        //
//                                                  //
//  Constructs the basis for the Sullivan space     //
//                                                  //
//  alpha is a multiindex                           //
//  sigma 1:k -> 0:n                                //
//  min[alpha] notin [sigma]                        //
//  [alpha] u [sigma] = [0..n]                      //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////




std::vector< std::pair<MultiIndex,IndexMap> > ListOfSullivanIndices( int n, int k, int r );








#endif
