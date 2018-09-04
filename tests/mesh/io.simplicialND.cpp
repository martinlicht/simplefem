

/**/

#include <iostream>
#include <sstream>
#include "../../basic.hpp"
#include "../../mesh/mesh.simplicialND.hpp"
#include "../../mesh/io.simplicialND.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for Simplicial ND Mesh IO" << endl;
    
    {
        
        MeshSimplicialND mesh = UnitSquareND();
        
        mesh.check();
        
        cout << "start IO..." << std::endl;
        
        std::stringstream ss;
        
        writeMeshSimplicialND( ss, mesh );
        
        std::cout << ss.str().c_str() << nl;
        
        ss.seekg( std::ios_base::beg );
        
        MeshSimplicialND mesh2 = readMeshSimplicialND( ss );
        
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
