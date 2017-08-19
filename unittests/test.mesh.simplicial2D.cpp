

/**/

#include <iostream>
#include <fstream>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/mesh.simplicial2D.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Simplicial 2D Module" << endl;
        
        MeshSimplicial2D M = TetrahedralSurface();
        
        M.check();
        
        cout << "Refinement" << endl;
        
        for( int c = 0; c < 2; c++ )
          M.uniformrefinement();
        
        M.check();
        
        cout << M << endl;
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
