#ifndef INCLUDEGUARD_GENERATEMULTIINDICES
#define INCLUDEGUARD_GENERATEMULTIINDICES

#include <algorithm>
#include <vector>
#include <iterator>

#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"
#include "multiindex.hpp"


/***************
*** 
***  Generate Multiindex lists 
***  
***************/

std::vector<MultiIndex> generateMultiIndices( const IndexRange&, int );



#endif