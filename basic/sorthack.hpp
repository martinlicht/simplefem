#ifndef INCLUDEGUARD_SORTHACK_HPP
#define INCLUDEGUARD_SORTHACK_HPP

#include <cassert>
#include <utility>
#include <vector>

template<typename A>
void sorthack( std::vector<A>& vec )
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
    
    for( int s = 0; s < vec.size(); s++ )
    for( int t = s; t > 0 && vec.at(t-1) > vec.at(t); t-- )
        std::swap( vec[t-1], vec[t] );
    
    for( int t = 1; t < vec.size(); t++ )
        assert( vec[t-1] <= vec[t] );
}

#endif
