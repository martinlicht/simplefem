
#include <iostream>
#include <new>

using namespace std;

/*
The purpose of this program is to check whether your development environment is set up correctly,
including compiler, linker, debug software, etc...
It also outputs the C++ version.
*/

int main( int argc, char *argv[] )
{
    cout << "Hello World! " << endl;
                
    cout << "C++ version: " << __cplusplus << endl;
        
    #ifdef __GLIBCXX__
    std::cout << "GNU libstdc++ version " << __GLIBCXX__ << std::endl;
    #endif

    #ifdef _LIBCPP_VERSION
    std::cout << "LLVM libc++ version " << _LIBCPP_VERSION << std::endl;
    #endif

    #ifdef _MSC_VER
    std::cout << "MSVC standard library with _MSC_VER=" << _MSC_VER << std::endl;
    #endif

    cout << "Now an intentional leak..." << endl;
    
    int * p = new (std::nothrow) int[10000];
    
    std::cout << p[8] << std::endl;
    
    return 0;
}

