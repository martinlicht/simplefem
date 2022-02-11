// g++ sorthacktest.cpp -o sorthacktest

#include <vector>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include "../../basic.hpp"



int main()
{
    
    LOG << "This is a line" << nl;
    NOTE "This is a note";
    PING;
    PING;
    PING;
    LOG << "Dies ist ein test!" << 5 << nl;
    LOGPRINTF( "%i%c%d\n", 1, '-', 3 );
    PING;
    // openContext();
    NOTE "";
    NOTE "Notice";
    WARNING "Warning";
    PING;
    // closeContext();
    return 0;
}
