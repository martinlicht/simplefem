

/**/

#include <iostream>
#include <fstream>

#include "../basic.hpp"
#include "coordinates.hpp"
#include "simplicialmesh.hpp"
#include "generatesimplicialmesh.hpp"
#include "vtkwriter.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Simplicial Mesh" << endl;
	
	{
		
		// SimplicialMesh M = UnitCubeTriangulation(3,3);
		SimplicialMesh M = UnitCubeTriangulation3D(3,3,3);
		
		cout << M << endl;

                cout << "Print VTK-type file" << endl;
                
                fstream fs( "./gitter.vtk", std::fstream::out );
                
                {
                    
                    VTK_MeshWriter vtk( M, fs );
                    vtk.writePreamble( "Mein erster Test" );
                    vtk.writeCoordinateBlock();
                    vtk.writeTopDimensionalCells();
                    
                }
                
                fs.close();

        }
    
        
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
