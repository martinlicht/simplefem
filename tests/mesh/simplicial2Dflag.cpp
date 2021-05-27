

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
        LOG << "Unit Test for Simplicial 2D Module" << endl;
        
//        MeshSimplicial2D M = TetrahedralSurface2D();
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        M.automatic_dirichlet_flags();

        M.check();
            
        M.check_dirichlet_flags();
            
        LOG << "Refinement" << endl;
        
        for( int c = 0; c < 2; c++ )
          M.uniformrefinement();
        
        M.check();
        
//         M.check_dirichlet_flags();
            
        LOG << M << endl;
        
        M.outputTikZ( std::cout );
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
