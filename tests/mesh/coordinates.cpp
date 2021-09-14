

/**/

#include <iostream>
#include <sstream>
#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/io.coordinates.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Coordinates" );

int main()
{
	LOG << "Unit Test: " << TestName << endl;
    // LOG << "Unit Test for Coordinates" << endl;
    
    {
        
        Coordinates coords(5,0);
        
        assert( coords == coords );
        
        std::stringstream ss;
        
        writeCoordinates( ss, coords );
        
        ss.seekg( std::ios_base::beg );
        
        Coordinates coords2 = readCoordinates( ss );
        
        coords.check();
        coords2.check();
        assert( coords == coords2 );
    }
    
    LOG << "Finished Unit Test: " << TestName << endl;
    
    return 0;
}
