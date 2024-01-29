
#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for one-dimensional simplicial mesh" << nl;

    MeshSimplicial1D M = StandardInterval1D();

    LOG << "Check" << nl;

    M.check();

    LOG << M << nl;

    LOG << "Start refinement" << nl;

    for( int c = 0; c < 100; c++ ) 
    {
        int e = c % M.count_edges();
        LOG << "Bisect edge: " << e << nl;
        M.bisect_edge( e );

        M.check();
    }

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
