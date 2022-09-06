/**/

#include <ostream>
// #include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicialND.hpp"
#include "../../mesh/examplesND.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "N-dimensional simplicial mesh" );

int main()
{
	LOG << "Unit Test: " << TestName << endl;
    
    WARNING "NOTHING IMPLEMENTED YET";
    
    if(false)
    {
        MeshSimplicialND M = HypertetrahedralSurface4D();
        
        LOG << "Check" << endl;
        
        M.check();
	
        LOG << "Check done" << endl;
        
        LOG << M << endl;
        
    }

	LOG << "Finished Unit Test: " << TestName << endl;
    
    return 0;
}
