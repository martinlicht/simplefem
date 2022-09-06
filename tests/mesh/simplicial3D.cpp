

/**/

#include <ostream>
// #include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Simplicial 3D Module" );

int main()
{
	LOG << "Unit Test: " << TestName << endl;
        // LOG << "Unit Test for Simplicial 3D Module" << endl;
        
        MeshSimplicial3D M = UnitSimplex3D();
        
        M.check();
        
        LOG << "Refinement..." << endl;
        
        M.uniformrefinement();
        
        LOG << "...done" << endl;
        
        M.check();
        
        LOG << M << endl;
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
