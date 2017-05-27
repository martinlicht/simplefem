

/**/

#include <iostream>
#include <fstream>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/manifold.2D.hpp"
#include "../mesh/vtkwriter.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Manifold 2D Module" << endl;
	
        ManifoldTriangulation2D M;
	
        cout << M << endl;
        
        cout << "Finished Unit Test" << endl;

	return 0;
}
