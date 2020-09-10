

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"

using namespace std;

int main()
{
        cout << "Unit Test for Simplicial 2D Module" << endl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        cout << "Refinement" << endl;
        
        for( int c = 0; c < 5; c++ ) {

          M.midpoint_refinement_global();
          
        }
          
        
        M.check();
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
