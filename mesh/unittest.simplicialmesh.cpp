

/**/

#include <iostream>
#include <fstream>

#include "../basic.hpp"
#include "coordinates.hpp"
#include "simplicialmesh.hpp"
#include "generatesimplicialmesh.hpp"
#include "vtkwriter.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Simplicial Mesh" << endl;
	
        SimplicialMesh M = UnitCubeTriangulation3D(3,3,3);
	
        cout << M << endl;
        
        cout << "Finished Unit Test" << endl;

	return 0;
}
