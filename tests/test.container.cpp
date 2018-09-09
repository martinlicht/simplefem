#include <iostream>

#include "../container.hpp"

using namespace std;

int main() {

    std::clog << "Unit test Container class!" << std::endl;
    
    ContainerInterface<double,int> foo(
        []() -> int { clog << "create iterator" << nl; return 0; },
        [](int&) -> void { clog << "destroy iterator" << nl;  },
        [](const int& i) -> bool { clog << "check end-condition" << nl; return i<10; },
        [](int& i) -> void { clog << "increment" << nl; i++; },
        [](const int& i) -> double { clog << "dereference" << nl; return i*i; }
    );
    
    for( auto it : foo )
      std::clog << it << nl;
    
    return 0;
}