

/**/

#include <iostream>
#include <fstream>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/mesh.simplicial2D.hpp"
#include "../vtk/vtkwriter.mesh3D.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for VTK output of Simplicial Mesh" << endl;
        
        {
                
                MeshSimplicial2D M = UnitSquare();
                // TODO: This should be a 3D mesh class
                
                cout << M << endl;
                
                cout << "Print VTK-type file" << endl;
                
                fstream fs( "./gitter.vtk", std::fstream::out );
                
                {
                    
                    VTK_MeshWriter_Mesh3D vtk( M, fs );
                    vtk.writePreamble( "Mein erster Test" );
                    vtk.writeCoordinateBlock();
                    vtk.writeTopDimensionalCells();
                    
                }
                
                fs.close();

        }
    
        
        
        cout << "Finished Unit Test" << endl;

        return 0;
}
