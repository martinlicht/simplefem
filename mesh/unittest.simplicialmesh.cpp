

/**/

#include <iostream>
#include "../basic.hpp"
#include "coordinates.hpp"
#include "simplicialmesh.hpp"
#include "generatesimplicialmesh.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Simplicial Mesh" << endl;
	
	{
		
		SimplicialMesh M = UnitCubeTriangulation(2,2);
		
		cout << M << endl;
		
	}
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
