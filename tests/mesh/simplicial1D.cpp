

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test for one-dimensional simplicial mesh" << endl;

        MeshSimplicial1D M = StandardInterval1D();
        
        LOG << "Check" << endl;
        
        M.check();
        
        LOG << M << endl;
        
        LOG << "Start refinement" << endl;
        
        for( int c = 0; c < 100; c++ ) 
        {
          int e = c % M.count_edges();
          LOG << "Bisect edge: " << e << nl;
          M.bisect_edge( e );
          
          M.check();
        }
        
        LOG << "Finished Unit Test" << endl;

        return 0;
}
