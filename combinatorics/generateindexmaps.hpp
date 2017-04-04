#ifndef INCLUDEGUARD_GENERATEINDEXMAPS
#define INCLUDEGUARD_GENERATEINDEXMAPS

#include <algorithm>
#include <vector>
#include <iterator>

#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"


/***************
*** 
***  Generate Index Maps of different kinds 
***  0) generate empty mapping 
***      1) generate all possible mappings 
***  2) generate Permuations, and tell their sign 
***  3) generate the Sigma mappings 
***  
***************/

std::vector<IndexMap> generateEmptyMap( const IndexRange& from, const IndexRange& to );

std::vector<IndexMap> generateIndexMaps( const IndexRange& from, const IndexRange& to );

std::vector<IndexMap> generatePermutations( const IndexRange& ir );

// std::vector<IndexMap> generateTransposition( const IndexRange& ir, int uno, int duo );

int signPermutation( const IndexMap& im );

std::vector<IndexMap> generateSigmas( const IndexRange& from, const IndexRange& to );



#endif