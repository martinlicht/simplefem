

/**/

#include <iostream>
#include <fstream>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/mesh.manifold1D.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Manifold 1D Module" << endl;
	
        MeshManifold1D M = UnitSquare();
        
        cout << "Check" << endl;
        
        M.check();
	
        cout << M << endl;
        
        cout << "Start refinement" << endl;
        
        M.improved_uniformrefinement();
        
        cout << "Finished Unit Test" << endl;

	return 0;
}
