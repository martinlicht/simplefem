

/**/

#include <ostream>
// #include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Manifold 1D Module" );

int main()
{
	LOG << "Unit Test: " << TestName << endl;
    // LOG << "Unit Test for Manifold 1D Module" << endl;

    MeshSimplicial1D M = StandardInterval1D();
    
    LOG << "Check" << endl;
    
    M.check();

    LOG << M << endl;
    
    LOG << "Start refinement" << endl;
    
    for( int c = 0; c < 10; c++ )
        M.improved_uniformrefinement();
    
    LOG << "Finished Unit Test: " << TestName << endl;

    return 0;
}
