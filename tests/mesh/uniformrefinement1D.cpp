

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Manifold 1D Module" << endl;
	
        MeshSimplicial1D M = UnitSquare1D();
        
        cout << "Check" << endl;
        
        M.check();
	
        cout << M << endl;
        
        cout << "Start refinement" << endl;
        
        for( int c = 0; c < 10; c++ )
          M.improved_uniformrefinement();
        
        cout << "Finished Unit Test" << endl;

	return 0;
}
