#ifndef INCLUDEGUARD_COMBINATORICS_HEAPS_ALGORITHM
#define INCLUDEGUARD_COMBINATORICS_HEAPS_ALGORITHM


#include <vector>

#include "../basic.hpp"


/***************
*** 
***  Generate Permutations using Heaps algorithm
***  
***  Init - initiliazes the auxiliary data: i and c
***  
***  Step - Returns true if a has experienced a transposition
***         Returns false if a has not been changed.
***         Supposed to be used in a do-while loop.
***  
***************/

void HeapsAlgorithmInit( int& i, std::vector<int>& c, const std::vector<int>& a );

bool HeapsAlgorithmStep( int& i, std::vector<int>& c, std::vector<int>& a );


#endif
