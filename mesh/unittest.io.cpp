

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
    
    cout << "Loading Unit Cube" << endl;
    SimplicialMesh mesh = UnitCubeTriangulation( 3, 3 );
    
    cout << "Writing and Reading" << endl;
    {

        cout << "Writing..." << endl;
        writeSimplicialMeshPath( "unittest.meshio.temp", mesh, false );
        
        cout << "Reading..." << endl;
        auto newmesh = readSimplicialMeshPath( "unittest.meshio.temp" );
        
        cout << "Writing on screen (1/2)..." << endl;
        writeSimplicialMesh( cout, mesh, true );
        
        cout << "Writing on screen (2/2)..." << endl;
        writeSimplicialMesh( cout, newmesh, true );
        
    }


    cout << "Finished Unit Test" << endl;

    return 0;
}
