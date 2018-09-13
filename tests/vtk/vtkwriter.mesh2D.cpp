

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../vtk/vtkwriter.mesh2D.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for VTK output of Simplicial Mesh" << endl;
    
    {
        
        // MeshSimplicial2D M = UnitCubeTriangulation(3,3);
        MeshSimplicial2D M = UnitSquare2D();
        
        cout << M << endl;

                cout << "Print VTK-type file" << endl;
                
                fstream fs( "./gitter.vtk", std::fstream::out );
                
                {
                    
                    VTK_MeshWriter_Mesh2D vtk( M, fs );
                    vtk.writePreamble( "Mein erster Test" );
                    vtk.writeCoordinateBlock();
                    vtk.writeTopDimensionalCells();
                    
                }
                
                fs.close();

        }
    
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
