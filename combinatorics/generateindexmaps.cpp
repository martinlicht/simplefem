
#include "generateindexmaps.hpp"

#include <algorithm>
#include <vector>
#include <iterator>

#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"


std::vector<IndexMap> 
generateEmptyMap( const IndexRange& from, const IndexRange& to )
{
    
    from.check();
    to.check();
    assert( from.isempty() );
    
    IndexMap myempty = IndexMap( from, to );
    myempty.check();
    std::vector<IndexMap> ret;
    ret.push_back( myempty );
    return ret;
    
}


std::vector<IndexMap> 
generateIndexMaps( const IndexRange& from, const IndexRange& to )
{
	
        from.check();
	to.check();
	
        if( from.isempty() )
            return generateEmptyMap( from, to );
        
        if( to.isempty() )
            return std::vector<IndexMap>();
        
        assert( from.cardinality() > 0 );
        assert( to.cardinality() > 0 );
        
        const int num = integerpower( to.cardinality(), from.cardinality() );
	std::vector<IndexMap> ret;
//         ( num, IndexMap( from, to ) );
	for( int it = 0; it < num; it++ ) {
                
                std::vector<int> values( from.getlength() );
                int it_clone = it;
                for( int digit = 0; digit < from.getlength(); digit++ ) {
                    values[digit] = to.min() + it_clone % to.cardinality();
                    it_clone /= to.cardinality();
                }
		/*
                for( int digit = 0; digit < from.getlength(); digit++ ) {
			ret[it][ from.min() + digit ] 
			=
			to.min() + 
			( it / integerpower( to.getlength(), digit ) ) % ( to.getlength() );
		}*/
		
		ret.push_back( IndexMap( from, to, values ) );
		
	}
	
	for( const auto& foo : ret )
            assert( (foo.check(),1) );
        
        return ret;
        
}

std::vector<IndexMap> 
generatePermutations( const IndexRange& ir )
{
	
        ir.check();
	std::vector<IndexMap> allmappings = generateIndexMaps( ir, ir );
	std::vector<IndexMap> ret( factorial( ir.getlength() ), IndexMap(ir) );
	std::copy_if( allmappings.begin(), allmappings.end(), ret.begin(),
		[]( const IndexMap& im ) -> bool { return im.isbijective(); }
		);
        
        for( const auto& perm : ret )
            assert( (perm.check(),1) && perm.isbijective() );
        
        return ret;
        
}

int signPermutation( const IndexMap& im )
{
	
    im.check();
    assert( im.isbijective() );
    assert( im.getSourceRange() == im.getDestRange() );
    
    const IndexRange& ir = im.getSourceRange();
    int zaehler = 1;
    int nenner = 1;
    
    for( int s = ir.min(); s <= ir.max(); s++ )
    for( int t = s+1; t <= ir.max(); t++ )
    {
        nenner *= ( t - s );
        zaehler *= ( im[ t ] - im[ s ] );
    }
    
    assert( zaehler % nenner == 0 );
    int ret = zaehler / nenner;
    assert( ret == 1 || ret == -1 );
    return ret;
    
}

std::vector<IndexMap> 
generateSigmas( const IndexRange& from, const IndexRange& to )
{
	from.check();
	to.check();
	std::vector<IndexMap> allmappings = generateIndexMaps( from, to );
	std::vector<IndexMap> ret( binomial<int>( to.getlength(), from.getlength() ), IndexMap(from,to) );
	std::copy_if( allmappings.begin(), allmappings.end(), ret.begin(),
		[]( const IndexMap& im ) -> bool { return im.isstrictlyascending(); }
		);
	return ret;
}



