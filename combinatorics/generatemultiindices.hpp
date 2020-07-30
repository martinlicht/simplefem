#ifndef INCLUDEGUARD_COMBINATORICS_GENERATEMULTIINDICES
#define INCLUDEGUARD_COMBINATORICS_GENERATEMULTIINDICES


#include <algorithm>
#include <iterator>
#include <vector>

#include "../basic.hpp"

#include "indexrange.hpp"
#include "indexmap.hpp"
#include "multiindex.hpp"


/***************
*** 
***  Generate Multiindex lists 
***  
***************/

std::vector<MultiIndex> generateMultiIndices( const IndexRange& ir, int absval );



#endif
