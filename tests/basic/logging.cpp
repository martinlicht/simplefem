// g++ sorthacktest.cpp -o sorthacktest

#include <vector>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "../../basic.hpp"



int main()
{
    
    PING;
    PING;
    PING;
    LOG << "Dies ist ein test!" << 5;
    PING;
    openContext();
    WARNING "Hallo Welt!";
    PING;
    NOTE "Hallo Welt!";
    closeContext();
    return 0;
}
