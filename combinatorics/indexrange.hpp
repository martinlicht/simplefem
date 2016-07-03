#ifndef INCLUDEGUARD_INDEXRANGE
#define INCLUDEGUARD_INDEXRANGE


#include <vector>
#include <limits>
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
		int getlength() const;
		
		bool isempty() const;
		bool contains( int ) const;
		bool contains( const IndexRange& subir ) const;
		bool operator== ( const IndexRange& ) const;
		int place( int ) const;
		
	private:

		int low;
		int high;
		
};

inline std::ostream& operator<<( std::ostream& os, const IndexRange& ir )
{
	ir.print( os );
	return os;
}

static const IndexRange  Nat0 = IndexRange( 0, std::numeric_limits<int>::max() );
static const IndexRange  Nat1 = IndexRange( 1, std::numeric_limits<int>::max() );
static const IndexRange& Nat  = Nat1;




#endif