#ifndef INCLUDEGUARD_GENERATEINDEXMAPS
#define INCLUDEGUARD_GENERATEINDEXMAPS

#include <algorithm>
#include <vector>
#include <iterator>

#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"


std::vector<IndexMap> generateIndexMaps( const IndexRange& from, const IndexRange& to );

std::vector<IndexMap> generatePermutations( const IndexRange& ir );

std::vector<IndexMap> generateSigmas( const IndexRange& from, const IndexRange& to );



#endif