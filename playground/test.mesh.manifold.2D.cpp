

/**/

#include <iostream>
#include <fstream>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/manifold.2D.hpp"


using namespace std;

int main()
{
	clog << "Unit Test for Manifold 2D Module" << endl;
	
        ManifoldTriangulation2D M = UnitSquare();
        
        M.check();
	
        clog << M << endl;
        
        clog << "Start refinement" << endl;
        
        for( int c = 0; c < 100; c++ ) 
        {
          int e = c % M.count_edges();
          clog << "Bisect edge: " << e << nl;
          M.bisect_edge( e );
          
//           clog << M << endl;
          
          M.check();
        }
        
        clog << "Finished Unit Test" << endl;

	return 0;
}
