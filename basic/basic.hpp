#ifndef INCLUDEGUARD_BASIC_HPP
#define INCLUDEGUARD_BASIC_HPP

#include <ctime>     
#include <cstdlib>     
#include <cassert>     /* assert macro */
#include <list>
#include <iterator>
#include <functional>
#include <iostream>

// // #include "assertion.hpp"

typedef double Float;

static const char space = ' ';

static const char* emptystring = "";

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
        { assert(false); return n; }
    else
        return n * factorial<T>(n-1);
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
T maximum( const T& a, const T& b )
{
    if( a >= b )
        return a;
    else
        return b;
}

template<typename T>
T minimum( const T& a, const T& b )
{
    if( a <= b )
        return a;
    else
        return b;
}

int binomial(int,int);


template<typename T>
T binomial( const T& n, const T& k )
{
    assert( 0 <= n );
    assert( 0 <= k && k <= n );
    return factorial(n) / ( factorial(k) * factorial(n-k) );
}


static inline int integerpower( int base, int exponent )
{
    assert( exponent != 0 || base != 0 );
    assert( exponent >= 0 );
    if( exponent == 0 ) return 1;
    return base * integerpower( base, exponent - 1 );
}

static inline int getbit( unsigned int value, unsigned int bitnumber )
{
    return ( value >> bitnumber ) % 2;
}





typedef clock_t timestamp;

inline timestamp gettimestamp()
{
    return CLOCKS_PER_SEC * static_cast<double>(clock()); 
}


template<typename T>
void mergeelementsinsortedlist
( std::list<T>& L, 
  std::function<T( const T&, const T& )> merge,
  std::function<bool( const T&, const T& )> compare
) {
    typename std::list<T>::iterator it = L.begin();
    while( it != L.end() ){

        typename std::list<T>::iterator now = it; 
        typename std::list<T>::iterator next = ++it;

        if( next == L.end() ) return;

        if( compare( *it, *next ) )
        {
            *now = merge( *now, *next );
            L.erase( next );
            it = now;
        } 

    }
}





inline void ping() { std::cout << "ping" << std::endl; }

inline void pong() { std::cout << "pong" << std::endl; }

inline void peng() { std::cout << "peng" << std::endl; }

inline void pang() { std::cout << "pang" << std::endl; }

inline void pung() { std::cout << "pung" << std::endl; }









#endif
