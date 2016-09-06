#ifndef INCLUDEGUARD_MULTIINDEX
#define INCLUDEGUARD_MULTIINDEX


#include <vector>
#include <iostream>
#include "../basic.hpp"
#include "indexrange.hpp"


/********************
**** 
****  Class that describes multiindices over an IndexRange 
****  basic arithmetic functionality
**** 
********************/

class MultiIndex
{

	public:
	
		MultiIndex( const IndexRange& ir );
		MultiIndex( const IndexRange& ir, const std::vector<int>& );
		
		void check() const;
		void print( std::ostream& ) const;
		
		IndexRange getIndexRange() const;
		const int& operator[](int) const;
		int& operator[](int);
		
		int absolute() const;
		int factorial() const;
		
		void operator+=(int);
		void operator+=( const MultiIndex&);
		void operator-=(int);
		void operator-=( const MultiIndex&);
		
		bool less( const MultiIndex& ) const;
		bool equals( const MultiIndex& ) const;
		
	private:

		IndexRange range;
		std::vector<int> values;
		
};




inline std::ostream& operator<<( std::ostream& os, const MultiIndex& mi )
{
	mi.check();
	mi.print( os );
	return os;
}

inline MultiIndex operator+( const MultiIndex& left, int right )
{
	left.check();
	MultiIndex ret = left;
	ret += right;
	ret.check();
	return ret;
}

inline MultiIndex operator-( const MultiIndex& left, int right )
{
	left.check();
	MultiIndex ret = left;
	ret -= right;
	ret.check();
	return ret;
}

inline MultiIndex operator+( const MultiIndex& left, const MultiIndex& right )
{
	left.check();
	right.check();
	MultiIndex ret = left;
	ret += right;
	ret.check();
	return ret;
}

inline MultiIndex operator-( const MultiIndex& left, const MultiIndex& right )
{
	left.check();
	right.check();
	MultiIndex ret = left;
	ret -= right;
	ret.check();
	return ret;
}


inline bool operator==( const MultiIndex& it, const MultiIndex& mi)
{
	it.check();
	mi.check();
	return it.equals( mi );
}
		
inline bool operator!=( const MultiIndex& it, const MultiIndex& mi)
{
	it.check();
	mi.check();
	return ! ( it == mi );
}

inline bool operator<( const MultiIndex& it, const MultiIndex& mi)
{
	it.check();
	mi.check();
	return it.less( mi );
}
		
inline bool operator>( const MultiIndex& it, const MultiIndex& mi)
{
	it.check();
	mi.check();
	return mi < it;
}
		
inline bool operator<=( const MultiIndex& it, const MultiIndex& mi)
{
	it.check();
	mi.check();
	return it < mi || it == mi;
}
		
inline bool operator>=( const MultiIndex& it, const MultiIndex& mi)
{
	it.check();
	mi.check();
	return it > mi || it == mi;
}
	
	

		
#endif