
#include "../../base/include.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"

// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: Simplicial 3D Module" << nl;
    
    MeshSimplicial3D M = UnitSimplex3D();
    
    M.check();

    M.automatic_dirichlet_flags();

    M.check_dirichlet_flags();
    
    LOG << "Refinement" << nl;
    
    for( int c = 0; c < 4; c++ ) {

        LOG << "Level " << c << " with " << M.count_tetrahedra() << " tetrahedra." << nl; 
        M.midpoint_refinement_global();
        
    }
        
    M.check();

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
