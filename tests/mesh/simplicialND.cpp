

/**/

#include <ostream>
// #include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicialND.hpp"
#include "../../mesh/examplesND.hpp"


using namespace std;

int main()
{
	LOG << "Unit Test for N-dimensional simplicial mesh" << endl;
	
    WARNING "NOTHING IMPLEMENTED YET";
    
    return 0; // TODO: complete the class for N-dimensional meshes in the first place

    MeshSimplicialND M = HypertetrahedralSurface4D();
        
        LOG << "Check" << endl;
        
        M.check();
	
        LOG << "Check done" << endl;
        
        LOG << M << endl;
        
        LOG << "Finished Unit Test" << endl;

	return 0;
}
