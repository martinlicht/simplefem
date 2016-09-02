

/**/

#include <iostream>
#include "../basic.hpp"
#include "coordinates.hpp"
#include "simplicialmesh.hpp"
#include "generatesimplicialmesh.hpp"
#include "io.coordinates.hpp"
#include "io.simplicialmesh.hpp"


using namespace std;

int main()
{
    cout << "Input/Output of Coordinates" << endl;

    SimplicialMesh mesh = UnitCubeTriangulation( 3, 3 );
    
    {

        writeSimplicialMeshPath( "unittest.meshio.temp", mesh, false );
        
        auto newmesh = readSimplicialMeshPath( "unittest.meshio.temp" );
        
        writeSimplicialMesh( cout, mesh, true );
        
        writeSimplicialMesh( cout, mesh, true );
        
    }


    cout << "Finished Unit Test" << endl;

    return 0;
}
