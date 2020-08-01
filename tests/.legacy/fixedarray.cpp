// g++ sorthacktest.cpp -o sorthacktest

#include <vector>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "../../basic.hpp"



int main()
{
    
        FixedArray<int> foo( 3 );
        for( auto i : foo )
            std::cout << i << ' ';
        std::cout << std::endl;
        
        FixedArray<int> goo( 3 );
        for( const auto i : goo )
            std::cout << i << ' ';
        std::cout << std::endl;
        
        const FixedArray<int> hoo( 3, 77);
        for( const auto i : hoo )
            std::cout << i << ' ';
        std::cout << std::endl;
        
        const FixedArray<int> ioo( 3, [](int i)->int{return i*i;});
        for( const auto i : ioo )
            std::cout << i << ' ';
        std::cout << std::endl;
        
        const FixedArray<int> joo( 3, {2,3,5} );
        for( const auto i : joo )
            std::cout << i << ' ';
        std::cout << std::endl;
        
}
