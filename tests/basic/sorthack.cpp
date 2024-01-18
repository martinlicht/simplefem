// g++ sorthacktest.cpp -o sorthacktest

#include <cstdlib>
#include <cassert>

#include <vector>
#include <algorithm>
#include "../../basic/sorthack.hpp"



int main()
{
    int N = 100;
    std::vector<int> foo( N );
    for( int i = 0; i < N; i++ ) foo[i] = rand() / 100000;

    sorthack( foo );
    
    for( int i = 1; i < N; i++ ) 
        assert( foo[i] >= foo[i-1] ); //std::cout << foo[i] << nl;
    
    return 0;
}
