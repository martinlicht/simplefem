

/**/

#include <ostream>
#include <sstream>
#include "../../basic.hpp"
#include "../../mesh/mesh.simplicialND.hpp"
#include "../../mesh/io.simplicialND.hpp"
#include "../../mesh/examplesND.hpp"


using namespace std;

int main()
{
    LOG << "Unit Test for Simplicial ND Mesh IO" << nl;
    
    WARNING "NOTHING IMPLEMENTED YET";
    
    if(false)
    {
        
        MeshSimplicialND mesh = UnitSquareND();
        
        mesh.check();
        
        LOG << "start IO..." << nl;
        
        std::stringstream ss;
        
        writeMeshSimplicialND( ss, mesh );
        
        LOG << ss.str().c_str() << nl;
        
        ss.seekg( std::ios_base::beg );
        
        MeshSimplicialND mesh2 = readMeshSimplicialND( ss );
        
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
