
#include <iostream>

using namespace std;

int main()
{
        cout << "Hello World! " << __cplusplus << endl;
        
#if __cplusplus < 201703L
        cout << "The version is " << __cplusplus << endl;
#endif 
        
        return 0;
}
