
#include "heappermgen.hpp"

#include <cassert>
#include <utility>
#include <vector>

#include "../basic.hpp"



void HeapsAlgorithmInit( int& i, std::vector<int>& c, const std::vector<int>& a )
{
    i = 0;
    assert( c.size() == a.size() );
    for( int j = 0; j < c.size(); j++ ) c[i] = 0;
}


bool HeapsAlgorithmStep( int& i, std::vector<int>& c, std::vector<int>& a )
{
    const int n = c.size();
    assert( c.size() == a.size() );
    assert( 0 <= i && i <= n );

    while( i < n ) {
        
        if( c[i] < i ) {
            
            if( i % 2 == 0 ) {
                std::swap( a[  0 ], a[i] );
            } else {
                std::swap( a[c[i]], a[i] );
            }
            
            c[i]++;
            i = 0;
            return true;
            
        } else {
            
            c[i] = 0;
            i++;
            
        }
        
    }

    return false;
}



