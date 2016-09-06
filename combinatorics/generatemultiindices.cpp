
#include "generateindexmaps.hpp"

#include <algorithm>
#include <vector>
#include <iterator>

#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"
#include "multiindex.hpp"


std::vector<MultiIndex> generateMultiIndices( const IndexRange& ir, int degree )
{
    assert( degree >= 0 );
    ir.check();
    std::vector<MultiIndex> ret;
    
    if( ir.isempty() ) {
      ret.push_back( MultiIndex(ir) );
      return ret;
    }
    
    int max_candidate = integerpower( degree+1, ir.getlength() );
    int min_index = ir.min();
    int max_index = ir.max();
    
    for( int candidate = degree; candidate < max_candidate; candidate++ ) {
        
        MultiIndex mi_candidate( ir );
        int candidate_copy = candidate;
        
        for( int i = min_index; i <= max_index; i++ ) {
            mi_candidate[ i ] = candidate_copy % ( degree+1 );
            candidate_copy /= degree+1;
        }
        
        mi_candidate.check();
        if( mi_candidate.absolute() == degree )
            ret.push_back( mi_candidate );
        
    }
    
    assert( ret.size() == binomial<int>( ir.getlength()-1 + degree, ir.getlength()-1 ) );
    return ret;
}




