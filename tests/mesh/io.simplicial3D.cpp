

/**/

#include <ostream>
#include <sstream>
#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/io.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

int main()
{
    
    LOG << "Unit Test for Simplicial 3D Mesh IO" << nl;
    
    {
        
        MeshSimplicial3D mesh = UnitCube3D();
        
        mesh.check();
        
        LOG << "start IO..." << std::endl;
        
        std::stringstream ss;
        
        writeMeshSimplicial3D( ss, mesh );
        
        ss.seekg( std::ios_base::beg );
        
        LOG << "create next mesh..." << std::endl;
        
        MeshSimplicial3D mesh2 = readMeshSimplicial3D( ss );
        
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
