
#include "generateindexmaps.hpp"

#include <algorithm>
#include <vector>
#include <iterator>

#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"


std::vector<IndexMap> 
generateIndexMaps( const IndexRange& from, const IndexRange& to )
{
	const int num = integerpower( to.getlength(), from.getlength() );
	std::vector<IndexMap> ret( num, IndexMap( from, to ) );
	for( int it = 0; it < num; it++ ) {
		for( int digit = 0; digit < from.getlength(); digit++ ) {
			ret[it][ from.getlow() + digit ] 
			=
			( num / integerpower( to.getlength(), digit ) ) % to.getlength();
		}
	}
	return ret;
}



std::vector<IndexMap> 
generatePermutations( const IndexRange& ir )
{
	std::vector<IndexMap> allmappings = generateIndexMaps( ir, ir );
	std::vector<IndexMap> ret( factorial( ir.getlength() ), IndexMap(ir,ir) );
	std::copy_if( allmappings.begin(), allmappings.end(), ret.begin(),
		[]( const IndexMap& im ) -> bool { return im.isbijective();}
		);
	return ret;
}


std::vector<IndexMap> 
generateSigmas( const IndexRange& from, const IndexRange& to )
{
	std::vector<IndexMap> allmappings = generateIndexMaps( from, to );
	std::vector<IndexMap> ret( binomial<int>( to.getlength(), from.getlength() ), IndexMap(from,to) );
	std::copy_if( allmappings.begin(), allmappings.end(), ret.begin(),
		[]( const IndexMap& im ) -> bool { return im.isstrictlyascending();}
		);
	return ret;
}



