

/**/

#include <iostream>
#include <fstream>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/simplicialmesh.hpp"
#include "../mesh/generatesimplicialmesh.hpp"
#include "../mesh/vtkwriter.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Simplicial Mesh" << endl;
	
        SimplicialMesh M = UnitCubeTriangulation3D(3,3,3);
	
        cout << M << endl;
        
        cout << "Finished Unit Test" << endl;

	return 0;
}
