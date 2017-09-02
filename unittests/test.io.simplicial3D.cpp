

/**/

#include <iostream>
#include <sstream>
#include "../basic.hpp"
#include "../mesh/mesh.simplicial3D.hpp"
#include "../mesh/io.simplicial3D.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for Simplicial 3D Mesh IO" << endl;
    
    {
        
        MeshSimplicial3D mesh = StandardSquare3D();
        
        mesh.check();
        
        cout << "start IO..." << std::endl;
        
        std::stringstream ss;
        
        writeMeshSimplicial3D( ss, mesh );
        
        ss.seekg( std::ios_base::beg );
        
        MeshSimplicial3D mesh2 = readMeshSimplicial3D( ss );
        
        cout << "check original mesh..." << std::endl;
        mesh.check();
        cout << "check replicated mesh..." << std::endl;
        mesh2.check();
        cout << "check mesh equivalence..." << std::endl;
        assert( mesh == mesh2 );
    }
    
    cout << "Finished Unit Test" << endl;
    
    return 0;
}
