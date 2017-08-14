

/**/

#include <iostream>
#include <fstream>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/mesh.simplicial1D.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for one-dimensional simplicial mesh" << endl;
	
        MeshSimplicial1D M = UnitSquare();
        
        cout << "Check" << endl;
        
        M.check();
	
        cout << M << endl;
        
        cout << "Start refinement" << endl;
        
        for( int c = 0; c < 100; c++ ) 
        {
          int e = c % M.count_edges();
          cout << "Bisect edge: " << e << nl;
          M.bisect_edge( e );
          
          M.check();
        }
        
//         cout << M << endl;
          
        cout << "Finished Unit Test" << endl;

	return 0;
}
