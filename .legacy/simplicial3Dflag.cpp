
#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Simplicial 3D Module" << nl;

    MeshSimplicial3D M = UnitSimplex3D();

    M.check();

    LOG << "Set automatic Dirichlet flags..." << nl;

    M.automatic_dirichlet_flags();

    M.check();

    M.check_dirichlet_flags();

    LOG << "Refinement..." << nl;

    M.uniformrefinement();

    M.check();

    M.check_dirichlet_flags();

    LOG << "Standard output..." << nl;

    LOG << M << nl;

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
