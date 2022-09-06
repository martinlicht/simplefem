// g++ sorthacktest.cpp -o sorthacktest

#include <vector>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include "../../basic/sorthack.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Sort hack" );

int main()
{
    std::clog << "Unit Test: " << TestName << endl;
    
    int N = 100;
    std::vector<int> foo( N );
    for( int i = 0; i < N; i++ ) foo[i] = rand() / 100000;

    sorthack( foo );
    
    for( int i = 1; i < N; i++ ) 
        assert( foo[i] >= foo[i-1] ); //std::cout << foo[i] << std::endl;
    
    std::clog << "Finished Unit Test: " << TestName << endl;
    
    return 0;
}
