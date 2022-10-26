
#include <iostream>
#include <new>

using namespace std;

/*
The purpose of this program is to check whether your development environment is set up correctly,
including compiler, linker, debug software, etc...
It also outputs the C++ version.
*/

int main()
{
    cout << "Hello World! " << __cplusplus << endl;
                
    cout << "C++ The version is " << __cplusplus << endl;
        
    cout << "Now an intentional leak..." << endl;
    
    int * p = new (std::nothrow) int[10000];
    
    std::cout << p[8] << std::endl;
    
    return 0;
}

