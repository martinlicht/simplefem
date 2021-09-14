

/**/

#include <iostream>
#include <sstream>
#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/io.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Simplicial 1D Mesh IO" );

int main()
{
	LOG << "Unit Test: " << TestName << endl;
    // LOG << "Unit Test for Simplicial 1D Mesh IO" << endl;
    
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
    
    LOG << "Finished Unit Test: " << TestName << endl;
    
    return 0;
}
