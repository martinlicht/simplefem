#ifndef INCLUDEGUARD_NATMAP
#define INCLUDEGUARD_NATMAP


#include <vector>
#include <iostream>
#include <cassert>
#include "../basic.hpp"



class IndexRange
{

	public:
	
		IndexRange( int, int );
		
		void check() const;
		void print( std::ostream& ) const;
		
		int getlow() const;
		int gethigh() const;
		
		bool contains( int ) const;
		bool contains( const IndexRange& subir ) const;
		bool operator== ( const IndexRange& ) const;
		
	private:

		int low;
		int high;
		
};

inline std::ostream& operator<<( std::ostream& os, const IndexRange& ir )
{
	ir.print( os );
	return os;
}




class IndexMap
{
	
	private:
	
		IndexRange src;
		IndexRange dest;
	
		std::vector<int> values;
	
	public:
	
		IndexMap( IndexRange, IndexRange );
		
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

		IndexMap inverse() const;
	
		IndexMap skip( int ) const;
	
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