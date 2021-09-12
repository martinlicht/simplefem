
#include <iostream>
#include "../../basic.hpp"

using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Various arithmetic functionality" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;

        // check factorial limits
        // check a few computation of factorials
        // check binomials
        // check powers 
        
        LOG << "    largest number of which the factorial can be computed" << std::endl;

        LOG << "    case signed char       : " << largest_factorial_base<       signed char>() << std::endl;
        LOG << "    case signed short      : " << largest_factorial_base<      signed short>() << std::endl;
        LOG << "    case signed int        : " << largest_factorial_base<        signed int>() << std::endl;
        LOG << "    case signed long       : " << largest_factorial_base<       signed long>() << std::endl;
        LOG << "    case signed long long  : " << largest_factorial_base<  signed long long>() << std::endl;
        
        LOG << "    case unsigned char     : " << largest_factorial_base<     unsigned char>() << std::endl;
        LOG << "    case unsigned short    : " << largest_factorial_base<    unsigned short>() << std::endl;
        LOG << "    case unsigned int      : " << largest_factorial_base<      unsigned int>() << std::endl;
        LOG << "    case unsigned long     : " << largest_factorial_base<     unsigned long>() << std::endl;
        LOG << "    case unsigned long long: " << largest_factorial_base<unsigned long long>() << std::endl;
        
        for( int n = 0; n <= 120; n++ )
        {
            auto largest_base = largest_factorial_base_AUX( n, 2 );
            LOG << largest_base << "! = " << factorial_integer(largest_base) << " <= " << n << std::endl;
        }
            
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
