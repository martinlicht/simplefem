

/**/

#include <ostream>
#include <sstream>
#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/io.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"


using namespace std;

int main()
{
    LOG << "Unit Test for Simplicial 1D Mesh IO" << nl;
    
    {
        
        MeshSimplicial1D mesh = StandardInterval1D();
        
        for( int c = 0; c < 14; c++ ) mesh.improved_uniformrefinement();
        
        mesh.check();
        
        LOG << "start IO..." << std::endl;
        
        std::stringstream ss;
        
        writeMeshSimplicial1D( ss, mesh );
        
        ss.seekg( std::ios_base::beg );
        
        MeshSimplicial1D mesh2 = readMeshSimplicial1D( ss );
        
        LOG << "check original mesh..." << std::endl;
        mesh.check();
        LOG << "check replicated mesh..." << std::endl;
        mesh2.check();
        LOG << "check mesh equivalence..." << std::endl;
        assert( mesh == mesh2 );
    }
    
    LOG << "Finished Unit Test" << nl;
    
    return 0;
}
