
#include <sstream>
#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/io.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main()
{
    LOG << "Unit Test for Simplicial 2D Mesh IO" << nl;
    
    {
        
        MeshSimplicial2D mesh = StandardSquare2D();
        
        for( int c = 0; c < 3; c++ ) mesh.uniformrefinement();
        
        mesh.check();
        
        LOG << "start IO..." << nl;
        
        std::stringstream ss;
        
        writeMeshSimplicial2D( ss, mesh );
        
        ss.seekg( std::ios_base::beg );
        
        MeshSimplicial2D mesh2 = readMeshSimplicial2D( ss );
        
        LOG << "check original mesh..." << nl;
        mesh.check();
        LOG << "check replicated mesh..." << nl;
        mesh2.check();
        LOG << "check mesh equivalence..." << nl;
        assert( mesh == mesh2 );
    }
    
    LOG << "Finished Unit Test" << nl;
    
    return 0;
}
