#ifndef INCLUDEGUARD_COMBINATORICS_GENERATEMULTIINDICES_HPP
#define INCLUDEGUARD_COMBINATORICS_GENERATEMULTIINDICES_HPP


#include <vector>

#include "../basic.hpp"

#include "indexrange.hpp"
#include "multiindex.hpp"




/***************
*** 
***  Generate Multiindex lists 
***  
***************/

std::vector<MultiIndex> generateMultiIndices( const IndexRange& ir, int absval );



#endif
