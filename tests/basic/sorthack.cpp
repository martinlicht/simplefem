// g++ sorthacktest.cpp -o sorthacktest

#include <cstdlib>
#include <cassert>

#include <vector>
#include <algorithm>
#include "../../basic/sorthack.hpp"

int main()
{
    std::size_t N = 2 << 22;
    std::vector<int> foo( N );
    for( std::size_t i = 0; i < N; i++ ) foo[i] = rand() / 100000;

    std::sort( foo.begin(), foo.end() ); 
    //sorthack( foo );
    
    for( std::size_t i = 1; i < N; i++ ) 
        assert( foo[i-1] <= foo[i] ); //std::cout << foo[i] << nl;
    
    return 0;
}
