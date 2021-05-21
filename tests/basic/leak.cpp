
#include <iostream>
#include <new>

using namespace std;

int main()
{
    int * p = new (std::nothrow) int[10000];
    
    std::cout << p[8] << std::endl;
    
    return 0;
}
