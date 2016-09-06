
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
    attest( from.isempty() );
    IndexMap empty = IndexMap( from, to );
	std::vector<IndexMap> ret( 1, empty );
	return ret;
}


std::vector<IndexMap> 
generateIndexMaps( const IndexRange& from, const IndexRange& to )
{
	const int num = integerpower( to.getlength(), from.getlength() );
	std::vector<IndexMap> ret( num, IndexMap( from, to ) );
	for( int it = 0; it < num; it++ ) {
		for( int digit = 0; digit < from.getlength(); digit++ ) {
			ret[it][ from.min() + digit ] 
			=
			to.min() + 
			( it / integerpower( to.getlength(), digit ) ) % ( to.getlength() );
		}
	}
	return ret;
}

std::vector<IndexMap> 
generatePermutations( const IndexRange& ir )
{
	std::vector<IndexMap> allmappings = generateIndexMaps( ir, ir );
	std::vector<IndexMap> ret( factorial( ir.getlength() ), IndexMap(ir) );
	std::copy_if( allmappings.begin(), allmappings.end(), ret.begin(),
		[]( const IndexMap& im ) -> bool { return im.isbijective(); }
		);
	return ret;
}

int signPermutation( const IndexMap& im )
{
	attest( im.isbijective() );
	attest( im.getSourceRange() == im.getDestRange() );
	const IndexRange& ir = im.getSourceRange();
	int zaehler = 1;
	int nenner = 1;
    
    for( int s = ir.min(); s <= ir.max(); s++ )
    for( int t = s+1; t <= ir.max(); t++ ) {
        nenner *= ( t - s );
        zaehler *= ( im[ t ] - im[ s ] );
    }
    attest( zaehler % nenner == 0 );
	int ret = zaehler / nenner;
	attest( ret == 1 || ret == -1 );
	return ret;
}

std::vector<IndexMap> 
generateSigmas( const IndexRange& from, const IndexRange& to )
{
	std::vector<IndexMap> allmappings = generateIndexMaps( from, to );
	std::vector<IndexMap> ret( binomial<int>( to.getlength(), from.getlength() ), IndexMap(from,to) );
	std::copy_if( allmappings.begin(), allmappings.end(), ret.begin(),
		[]( const IndexMap& im ) -> bool { return im.isstrictlyascending(); }
		);
	return ret;
}



