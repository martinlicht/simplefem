

#include "generatemultiindices.hpp"

#include <vector>

#include "../basic.hpp"

#include "indexrange.hpp"
#include "multiindex.hpp"
#include "generateindexmaps.hpp"



std::vector<MultiIndex> generateMultiIndices( const IndexRange& ir, int degree )
{
    assert( degree >= 0 );
    ir.check();
    std::vector<MultiIndex> ret;
    
    /* if the index range is empty, then only the empty mapping is returned */
    
    if( ir.isempty() ) {
      ret.push_back( MultiIndex(ir) );
      return ret;
    }
    
    /* otherwise, we iterate through a large list of candidates */
    /* and pick those that have the right properties            */
    
    ret.reserve( binomial_integer( ir.cardinality()-1 + degree, ir.cardinality()-1 ) );
    
    int max_candidate = power_integer( degree+1, ir.cardinality() );
    int min_index = ir.min();
    int max_index = ir.max();
    
    for( int candidate = degree; candidate < max_candidate; candidate++ ) {
        
        MultiIndex mi_candidate( ir );
        int candidate_copy = candidate;
        
        for( int i = min_index; i <= max_index; i++ ) {
            mi_candidate.at( i ) = candidate_copy % ( degree+1 );
            candidate_copy /= degree+1;
        }
        
        mi_candidate.check();
        
        if( mi_candidate.absolute() == degree )
            ret.push_back( mi_candidate );
        
    }
    
    assert( ret.size() == binomial_integer( ir.cardinality()-1 + degree, ir.cardinality()-1 ) );
    
    assert( ret.size() == ret.capacity() );
    
    return ret;
}




