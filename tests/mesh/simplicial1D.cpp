
#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/io.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for one-dimensional simplicial mesh" << nl;

    MeshSimplicial1D M = StandardInterval1D();

    M.check();

    M.automatic_dirichlet_flags();

    M.check_dirichlet_flags();

    LOG << M << nl;

    LOG << "Start refinement" << nl;

    for( int c = 0; c < 0; c++ ) 
        M.improved_uniformrefinement();
        
    if(false)for( int c = 0; c < 100; c++ ) 
    {
        int e = c % M.count_edges();
        LOG << "Bisect edge: " << e << nl;
        M.bisect_edge( e );

        M.check();
        M.check_dirichlet_flags();
    }

    {
        LOG << "start IO..." << nl;
        
        std::stringstream ss;
        
        writeMeshSimplicial1D( ss, M );
        
        ss.seekg( std::ios_base::beg );
        
        MeshSimplicial1D M2 = readMeshSimplicial1D( ss );
        
        LOG << "check mesh equivalence..." << nl;
        M2.check();
        assert( M.getcoordinates() == M2.getcoordinates() );
        assert( M == M2 );
    }

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
