

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
// #include "../../mesh/mesh.simplicialND.hpp"
#include "../../vtk/vtkwriter.mesh3D.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for VTK output of Simplicial Mesh" << endl;
        
        {
                
                MeshSimplicial3D M = StandardCube3D();
                // TODO: This should be a 3D mesh class
                
                cout << M << endl;
                
                cout << "Print VTK-type file" << endl;
                
                fstream fs( "./gitter.vtk", std::fstream::out );
                
                {
                    
                    VTK_MeshWriter_Mesh3D vtk( M, fs, "Mein erster Test" );
                    vtk.writeCoordinateBlock();
                    vtk.writeTopDimensionalCells();
                    
                }
                
                fs.close();

        }
    
        
        
        cout << "Finished Unit Test" << endl;

        return 0;
}
