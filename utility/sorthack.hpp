#ifndef INCLUDEGUARD_SORTHACK_HPP
#define INCLUDEGUARD_SORTHACK_HPP

#include <algorithm>
#include <utility>
#include <vector>

#include "../base/include.hpp"

// ============================================================================
// Sort functions
// Originally, this function was introduced to mitigate a defect in C++14
// that prevented std::sort to apply a std::vector of std::array<int,N>.
// ============================================================================

template<typename A>
inline void sorthack( std::vector<A>& vec )
{

    std::sort( vec.begin(), vec.end() );
    
    // quickSort( vec );
    
    // for( int s = 0; s < vec.size(); s++ )
    // for( int t = 0; t < vec.size(); t++ )
    // {
    //     if( s < t && vec[s] > vec[t] )
    //     {
    //         A temp = vec[t];
    //         vec[t] = vec[s];
    //         vec[s] = temp;
    //     }
    // }

    // // insertion sort
    // for( int s = 0; s < vec.size(); s++ )
    // for( int t = s; t > 0 && vec.at(t-1) > vec.at(t); t-- )
    // {
    //     A temp = vec[t-1];
    //     vec[t-1] = vec[t];
    //     vec[t] = temp;
    // }

    // for( int t = 1; t < vec.size(); t++ )
    //     assert( vec[t-1] <= vec[t] );
}



#endif
