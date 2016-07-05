#ifndef INCLUDEGUARD_INDEXMAP
#define INCLUDEGUARD_INDEXMAP


#include <cassert>

#include <vector>
#include <limits>
#include <iostream>

#include "../basic.hpp"
#include "indexrange.hpp"






class IndexMap
{
	
	private:
	
		IndexRange src;
		IndexRange dest;
	
		std::vector<int> values;
	
	public:
	
		IndexMap( IndexRange );
		IndexMap( IndexRange, IndexRange );
		IndexMap( IndexRange, std::vector<int> );
		IndexMap( IndexRange, IndexRange, std::vector<int> );
		
		void check() const;
		void print( std::ostream& ) const;
		
		const IndexRange& getSourceRange() const;
		const IndexRange& getDestRange() const;
	
		bool isinjective() const;
		bool issurjective() const;
		bool isbijective() const;
		bool isstrictlyascending() const;
                
		int& operator[]( int i );
		const int& operator[]( int i ) const;
		const std::vector<int>& getvalues() const;
		
		IndexMap inverse() const;
	
		IndexMap skip( int ) const;
		IndexMap attachbefore( int ) const;
	
};



inline IndexMap operator*( const IndexMap& leave, const IndexMap& enter )
{
	
	IndexRange src  = enter.getSourceRange();
	IndexRange dest = leave.getDestRange();
	
	assert( enter.getDestRange() == leave.getSourceRange() );
	
	IndexMap ret( src, dest );
	
	for( int i = src.getlow(); i <= src.gethigh(); i++ )
		ret[i] = leave[ enter[i] ];
		
	return ret;
	
}

inline std::ostream& operator<<( std::ostream& os, IndexMap im )
{
	im.print( os );
	return os;
}



static IndexMap identityIndexMap( int low, int high )
{
	IndexRange ir( low, high );
	IndexMap im( ir, ir );
	for( int i = low; i <= high; i++ )
		im[i] = i;
	return im;
}



#endif