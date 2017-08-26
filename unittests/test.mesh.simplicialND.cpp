

/**/

#include <iostream>
#include <fstream>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/mesh.simplicialND.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for N-dimensional simplicial mesh" << endl;
	
        MeshSimplicialND M = HypertetrahedralSurface4D();
        
        cout << "Check" << endl;
        
        M.check();
	
        cout << "Check done" << endl;
        
        cout << M << endl;
        
        cout << "Finished Unit Test" << endl;

	return 0;
}
