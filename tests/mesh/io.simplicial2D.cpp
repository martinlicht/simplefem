

/**/

#include <ostream>
#include <sstream>
#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/io.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main()
{
    LOG << "Unit Test for Simplicial 2D Mesh IO" << endl;
    
    {
        
        MeshSimplicial2D mesh = StandardSquare2D();
        
        for( int c = 0; c < 3; c++ ) mesh.uniformrefinement();
        
        mesh.check();
        
        LOG << "start IO..." << std::endl;
        
        std::stringstream ss;
        
        writeMeshSimplicial2D( ss, mesh );
        
        ss.seekg( std::ios_base::beg );
        
        MeshSimplicial2D mesh2 = readMeshSimplicial2D( ss );
        
        LOG << "check original mesh..." << std::endl;
        mesh.check();
        LOG << "check replicated mesh..." << std::endl;
        mesh2.check();
        LOG << "check mesh equivalence..." << std::endl;
        assert( mesh == mesh2 );
    }
    
    LOG << "Finished Unit Test" << endl;
    
    return 0;
}
