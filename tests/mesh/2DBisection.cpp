

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Simplicial 2D Module" << endl;
        
        MeshSimplicial2D M = UnitSquare2D();
        
        M.check();
        
        for( int c = 0; c < 5; c++ ) {
          
            M.bisect_edge( c % M.count_edges() );
          
        }
        
        M.check();
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
