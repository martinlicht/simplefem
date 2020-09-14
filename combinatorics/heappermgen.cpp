
#include "heappermgen.hpp"

#include <cassert>
#include <utility>
#include <vector>

#include "../basic.hpp"



void HeapsAlgorithmInit( int& i, std::vector<int>& memo, const std::vector<int>& perm )
{
    i = 0;
    assert( memo.size() == perm.size() );
    for( int j = 0; j < memo.size(); j++ ) memo[i] = 0;
}


bool HeapsAlgorithmStep( int& i, std::vector<int>& memo, std::vector<int>& perm )
{
    const int n = memo.size();
    assert( memo.size() == perm.size() );
    assert( 0 <= i && i <= n );

    while( i < n ) {
        
        if( memo[i] < i ) {
            
            if( i % 2 == 0 ) {
                std::swap( perm[     0 ], perm[i] );
            } else {
                std::swap( perm[memo[i]], perm[i] );
            }
            
            memo[i]++;
            i = 0;
            return true;
            
        } else {
            
            memo[i] = 0;
            i++;
            
        }
        
    }

    return false;
}



