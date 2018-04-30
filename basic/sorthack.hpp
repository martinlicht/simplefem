#ifndef INCLUDEGUARD_SORTHACK_HPP
#define INCLUDEGUARD_SORTHACK_HPP

#include <vector>
#include <cassert>

template<typename A>
void sorthack( std::vector<A>& vec )
{
    for( int t = 0; t < vec.size(); t++ )
    for( int s = 0; s < vec.size(); s++ )
    {
        if( t < s && vec[t] > vec[s] )
        {
            A temp = vec[t];
            vec[t] = vec[s];
            vec[s] = temp;
        }
    }
    
    for( int t = 1; t < vec.size(); t++ )
        assert( vec[t-1] < vec[t] );
}

#endif
