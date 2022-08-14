

/**/

#include <ostream>
#include <sstream>
#include "../../basic.hpp"
#include "../../mesh/mesh.simplicialND.hpp"
#include "../../mesh/io.simplicialND.hpp"
#include "../../mesh/examplesND.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Simplicial ND Mesh IO" );

int main()
{
	LOG << "Unit Test: " << TestName << endl;
    // LOG << "Unit Test for Simplicial ND Mesh IO" << endl;
    
    WARNING "NOTHING IMPLEMENTED YET";
    
    if(false)
    {
        
        MeshSimplicialND mesh = UnitSquareND();
        
        mesh.check();
        
        LOG << "start IO..." << std::endl;
        
        std::stringstream ss;
        
        writeMeshSimplicialND( ss, mesh );
        
        LOG << ss.str().c_str() << nl;
        
        ss.seekg( std::ios_base::beg );
        
        MeshSimplicialND mesh2 = readMeshSimplicialND( ss );
        
        LOG << "check original mesh..." << std::endl;
        mesh.check();
        LOG << "check replicated mesh..." << std::endl;
        mesh2.check();
        LOG << "check mesh equivalence..." << std::endl;
        assert( mesh == mesh2 );
    }
    
    LOG << "Finished Unit Test: " << TestName << endl;
    
    return 0;
}
