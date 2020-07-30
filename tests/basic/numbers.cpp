
#include <iostream>
#include "../../basic.hpp"

using namespace std;

int main()
{
        cout << "Unit test for a few arithmetic functions" << endl;
        
        // check factorial limits
        // check a few computation of factorials
        // check binomials
        // check powers 
        
        std::cout << "    sizes of classical signed integer types" << std::endl;
        std::cout << "    case signed char       : " << largest_factorial_base<       signed char>() << std::endl;
        std::cout << "    case signed short      : " << largest_factorial_base<      signed short>() << std::endl;
        std::cout << "    case signed int        : " << largest_factorial_base<        signed int>() << std::endl;
        std::cout << "    case signed long       : " << largest_factorial_base<       signed long>() << std::endl;
        std::cout << "    case signed long long  : " << largest_factorial_base<  signed long long>() << std::endl;
        
        std::cout << "    sizes of classical signed integer types" << std::endl;
        std::cout << "    case unsigned char     : " << largest_factorial_base<     unsigned char>() << std::endl;
        std::cout << "    case unsigned short    : " << largest_factorial_base<    unsigned short>() << std::endl;
        std::cout << "    case unsigned int      : " << largest_factorial_base<      unsigned int>() << std::endl;
        std::cout << "    case unsigned long     : " << largest_factorial_base<     unsigned long>() << std::endl;
        std::cout << "    case unsigned long long: " << largest_factorial_base<unsigned long long>() << std::endl;
        
        for( int k = 0; k < 100; k++ )
            std::cout << k << " " << largest_factorial_base_AUX( k, 2 ) << std::endl;
        
        return 0;
}
