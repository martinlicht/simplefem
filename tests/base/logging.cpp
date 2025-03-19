// g++ sorthacktest.cpp -o sorthacktest


#include <cstdlib>

#include "../../base/include.hpp"

// ============================================================================
// Tests the fundamental features of the logging framework.
// ============================================================================


int main( int argc, char *argv[] )
{
    LOG << "Unit Test: Logging" << nl;
    
    LOG << "This is a line" << nl;
    NOTE "This is a note";
    PING;
    PING;
    PING;
    LOG << "This is a test! " << 5 << nl;
    LOGPRINTF( "%i%c%d\n", 1, '-', 3 );
    LOGPRINTF( "double      0.123456789e-7: %e\n",   (double)     0.1234567890123456789e-7 );
    LOGPRINTF( "long double 0.123456789e-7: % Le\n", (long double)0.1234567890123456789e-7 );
    PING;
    NOTE "";
    NOTE "Notice";
    WARNING "Warning";
    PING;
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
