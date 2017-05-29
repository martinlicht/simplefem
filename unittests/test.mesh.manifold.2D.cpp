

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
	
        ManifoldTriangulation2D M = UnitSquare();
        
        M.check();
	
        cout << M << endl;
        
        cout << "Start refinement" << endl;
        
        for( int c = 0; c < 100; c++ ) 
        {
          int e = c % M.count_edges();
          cout << "Bisect edge: " << e << nl;
          M.bisect_edge( e );
          
//           cout << M << endl;
          
          M.check();
        }
        
        cout << "Finished Unit Test" << endl;

	return 0;
}
