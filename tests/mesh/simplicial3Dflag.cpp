
#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

int main()
{
    LOG << "Unit Test for Simplicial 3D Module" << nl;

    MeshSimplicial3D M = UnitSimplex3D();

    M.check();

    M.automatic_dirichlet_flags();

    M.check_dirichlet_flags();

    LOG << "Refinement..." << nl;

    M.uniformrefinement();

    LOG << "...done" << nl;

    M.check();

    M.check_dirichlet_flags();

    LOG << M << nl;

    LOG << "Finished Unit Test" << nl;

    return 0;
}
