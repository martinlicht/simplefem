#ifndef INCLUDEGUARD_COMBINATORICS_GENERATEMULTIINDICES_HPP
#define INCLUDEGUARD_COMBINATORICS_GENERATEMULTIINDICES_HPP


#include <vector>

#include "../base/include.hpp"

#include "indexrange.hpp"
#include "multiindex.hpp"




/***************
*** 
***  Generate Multiindex lists 
***  
***  Generate all the multi-indices over an index range ir
***  with all the absolute values being `degree`.
***  
***  
***************/

std::vector<MultiIndex> generateMultiIndices( const IndexRange& ir, int degree );



#endif
