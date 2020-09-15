// g++ sorthacktest.cpp -o sorthacktest

#include <vector>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "../../basic/sorthack.hpp"



int main()
{
    int N = 100;
    std::vector<int> foo( N );
    for( int i = 0; i < N; i++ ) foo[i] = rand() / 100000;
    
    sorthack( foo );
    
    for( int i = 1; i < N; i++ ) 
        assert( foo[i] >= foo[i-1] ); //std::cout << foo[i] << std::endl;
    
    return 0;
}
