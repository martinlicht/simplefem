#ifndef INCLUDEDGARD_BASIC_HPP
#define INCLUDEDGARD_BASIC_HPP

#include <ctime>     
#include <cstdlib>     
#include <cassert>     /* assert macro */



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


template<typename T>
T binomial( const T& n, const T& k )
{
	assert( 0 <= n );
	assert( 0 <= k && k <= n );
	return factorial(n) / ( factorial(k) * factorial(n-k) );
}


static int integerpower( int base, int exponent )
{
	assert( exponent >= 0 );
	if( exponent == 0 ) return 1;
	return base * integerpower( base, exponent - 1 );
}




typedef clock_t timestamp;

inline timestamp gettimestamp()
{
	return CLOCKS_PER_SEC * static_cast<double>(clock()); 
}


#endif
