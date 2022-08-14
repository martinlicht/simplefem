

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
        
        MeshSimplicial3D M = UnitSimplex3D();
        
        M.check();
        
        M.automatic_dirichlet_flags();

        M.check_dirichlet_flags();
        
        LOG << "Refinement..." << endl;
        
        M.uniformrefinement();
        
        LOG << "...done" << endl;
        
        M.check();
        
        M.check_dirichlet_flags();
        
        LOG << M << endl;
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
