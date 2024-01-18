#ifndef INCLUDEGUARD_SORTHACK_HPP
#define INCLUDEGUARD_SORTHACK_HPP

#include <cassert>
#include <utility>
#include <vector>

template<typename A>
inline void sorthack( std::vector<A>& vec )
{
//     for( int s = 0; s < vec.size(); s++ )
//     for( int t = 0; t < vec.size(); t++ )
//     {
//         if( s < t && vec[s] > vec[t] )
//         {
//             A temp = vec[t];
//             vec[t] = vec[s];
//             vec[s] = temp;
//         }
//     }
    
    /* insertion sort */
    for( int s = 0; s < vec.size(); s++ )
    for( int t = s; t > 0 && vec.at(t-1) > vec.at(t); t-- )
    // {
    //     A temp = vec[t-1];
    //     vec[t-1] = vec[t];
    //     vec[t] = temp;
    // }
        std::swap( vec[t-1], vec[t] );
    // TODO: understand why swap for arrays is not found here
    
    for( int t = 1; t < vec.size(); t++ )
        assert( vec[t-1] <= vec[t] );
}

#endif
