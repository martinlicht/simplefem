#include<iostream>
#include<math.h>
#include<unistd.h>
#include<iomanip>
using namespace std;
int main()
{
    std::cout << std::endl;
    std::cout << std::scientific;
    std::cout << std::setprecision(5);
    const int cmax = 100;
    for( int c = 0; c < cmax; c++ ) 
    {
    
        std::cout << "\r" << c << "/" << cmax << " " << 2*c*sqrt(c);
        std::cout.flush();
        sleep(1);
    }
    std::cout << std::endl;
}
