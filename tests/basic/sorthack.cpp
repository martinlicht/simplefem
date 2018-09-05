// g++ sorthacktest.cpp -o sorthacktest

#include <vector>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "sorthack.hpp"



int main()
{
    int N = 100;
    std::vector<int> foo( N );
    for( int i = 0; i < N; i++ ) foo[i] = rand() / 100000;
    
//     std::sort( foo.begin(), foo.end() );
    sorthack( foo );
    
    for( int i = 0; i < N; i++ ) std::cout << foo[i] << std::endl;
    
    return 0;
}