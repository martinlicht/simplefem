#include <iostream>

#include "../container.hpp"

using namespace std;

int main() {

    std::cout << "Unit test Container class!" << std::endl;
    
    ContainerInterface<double,int> foo(
        []() -> int { cout << "create iterator" << nl; return 0; },
        [](int&) -> void { cout << "destroy iterator" << nl;  },
        [](const int& i) -> bool { cout << "check end-condition" << nl; return i<10; },
        [](int& i) -> void { cout << "increment" << nl; i++; },
        [](const int& i) -> double { cout << "dereference" << nl; return i*i; }
    );
    
    for( auto it : foo )
      std::cout << it << nl;
    
    return 0;
}