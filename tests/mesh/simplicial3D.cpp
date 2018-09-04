

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Simplicial 3D Module" << endl;
        
        MeshSimplicial3D M = UnitSimplex3D();
        
        M.check();
        
        cout << "Refinement..." << endl;
        
        M.uniformrefinement();
        
        cout << "...done" << endl;
        
        M.check();
        
        cout << M << endl;
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
