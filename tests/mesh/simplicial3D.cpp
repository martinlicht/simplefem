

/**/

#include <ostream>
// #include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test for Simplicial 3D Module" << endl;
        
        MeshSimplicial3D M = UnitSimplex3D();
        
        M.check();
        
        LOG << "Refinement..." << endl;
        
        M.uniformrefinement();
        
        LOG << "...done" << endl;
        
        M.check();
        
        LOG << M << endl;
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
