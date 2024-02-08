
#include <sstream>
#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/io.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    
    LOG << "Unit Test for Simplicial 3D Mesh IO" << nl;
    
    {
        
        MeshSimplicial3D mesh = UnitCube3D();
        
        mesh.check();
        
        LOG << "start IO..." << nl;
        
        std::stringstream ss;
        
        writeMeshSimplicial3D( ss, mesh );
        
        ss.seekg( std::ios_base::beg );
        
        LOG << "create next mesh..." << nl;
        
        MeshSimplicial3D mesh2 = readMeshSimplicial3D( ss );
        
        LOG << "check original mesh..." << nl;
        mesh.check();
        LOG << "check replicated mesh..." << nl;
        mesh2.check();
        LOG << "check mesh equivalence..." << nl;
        assert( mesh == mesh2 );
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
