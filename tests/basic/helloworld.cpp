
#include <iostream>
#include <new>

using namespace std;

int main()
{
    cout << "Hello World! " << __cplusplus << endl;
                
    #if __cplusplus < 201703L
    cout << "The version is " << __cplusplus << endl;
    #endif 
        
    cout << "Now an intentional leak..." << endl;
    
    int * p = new (std::nothrow) int[10000];
    
    std::cout << p[8] << std::endl;
    
    return 0;
}

