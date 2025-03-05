
#include "heappermgen.hpp"

#include <utility>
#include <vector>

#include "../basic.hpp"



void HeapsAlgorithmInit( int& iter, std::vector<int>& memo, const std::vector<int>& perm )
{
    iter = 0;
    assert( memo.size() == perm.size() );
    for( int j = 0; j < memo.size(); j++ ) memo[iter] = 0;
}


bool HeapsAlgorithmStep( int& iter, std::vector<int>& memo, std::vector<int>& perm )
{
    const std::vector<int>::size_type n = memo.size();
    assert( memo.size() == perm.size() );
    assert( 0 <= iter && iter <= n );

    while( iter < n ) {
        
        if( memo[iter] < iter ) {
            
            if( iter % 2 == 0 ) {
                std::swap( perm[        0 ], perm[iter] );
            } else {
                std::swap( perm[memo[iter]], perm[iter] );
            }
            
            memo[iter]++;
            iter = 0;
            return true;
            
        } else {
            
            memo[iter] = 0;
            iter++;
            
        }
        
    }

    return false;
}



