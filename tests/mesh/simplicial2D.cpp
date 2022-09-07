

/**/

#include <iostream>
// #include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test for Simplicial 2D Module" << nl;
        
//        MeshSimplicial2D M = TetrahedralSurface2D();
        MeshSimplicial2D M = StandardSquare2D_strange14();
        
        M.check();
        
        LOG << "Refinement" << nl;
        
        for( int c = 0; c < 2; c++ )
          M.uniformrefinement();
        
        M.check();
        
        LOG << M << nl;
        
        std::cout << M.outputTikZ();
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}
