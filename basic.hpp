#ifndef INCLUDEDGARD_BASIC_HPP
#define INCLUDEDGARD_BASIC_HPP

#include <cassert>     /* assert macro */
// #include <cstdlib>     /* abs */


typedef double Float;

static const char space = ' ';

static const char nl = '\n';

static const char tab = '\t';

inline int kronecker( int i, int j )
{
    if( i == j )
        return 1;
    else
        return 0;
}


template<typename T>
T factorial( const T& n )
{
	if( n == 0 )
		return 1;
	else if( n < 0 )
		assert(false);
	else
		return factorial<T>(n-1);
}

template<typename T>
T absolute( const T& n )
{
	if( n >= 0 )
		return n;
	else
		return -n;
}

int binomial(int,int);



#endif
