#ifndef INCLUDEGUARD_MULTIINDEX
#define INCLUDEGUARD_MULTIINDEX


#include <vector>
#include <iostream>
#include <cassert>
#include "../basic.hpp"
#include "indexrange.hpp"



class MultiIndex
{

	public:
	
		MultiIndex( IndexRange ir );
		
		void check() const;
		void print( std::ostream& ) const;
		
		IndexRange getIndexRange() const;
		const int operator[](int) const;
		int operator[](int);
		
		int absolute() const;
		int factorial() const;
		
		void operator+=(int);
		void operator+=( const MultiIndex&);
		void operator-=(int);
		void operator-=( const MultiIndex&);
		
	private:

		IndexRange range;
		std::vector<int> values;
		
};




inline std::ostream& operator<<( std::ostream& os, const MultiIndex& mi )
{
	mi.print( os );
	return os;
}

MultiIndex operator+( const MultiIndex& left, int right )
{
	MultiIndex ret = left;
	ret += right;
	return ret;
}

MultiIndex operator-( const MultiIndex& left, int right )
{
	MultiIndex ret = left;
	ret -= right;
	return ret;
}

MultiIndex operator+( const MultiIndex& left, const MultiIndex& right )
{
	MultiIndex ret = left;
	ret += right;
	return ret;
}

MultiIndex operator-( const MultiIndex& left, const MultiIndex& right )
{
	MultiIndex ret = left;
	ret -= right;
	return ret;
}

		


#endif