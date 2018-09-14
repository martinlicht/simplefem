

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../vtk/vtkwriter.mesh1D.hpp"
#include "../../mesh/examples1D.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for VTK output of Simplicial Mesh" << endl;
    
    {
        
        // MeshSimplicial1D M = UnitCubeTriangulation(3,3);
        MeshSimplicial1D M = UnitSquare1D();
        
        cout << M << endl;
        
        cout << "Print VTK-type file" << endl;
        
        fstream fs( "./gitter.vtk", std::fstream::out );
        
        {
            
            VTK_MeshWriter_Mesh1D vtk( M, fs );
            vtk.writePreamble( "Mein erster Test" );
            vtk.writeCoordinateBlock();
            vtk.writeTopDimensionalCells();
            
        }
        
        fs.close();
        
    }
    
    cout << "Finished Unit Test" << endl;
    
    return 0;
}
