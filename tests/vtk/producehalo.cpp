

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../vtk/vtkwriter.mesh2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for VTK output of Simplicial Mesh" << endl;
    
    const int K = 4;
    const int L = 12;
    
    {
        
        MeshSimplicial2D M = Halo(K,L);
        
        cout << "Print VTK-type file" << endl;
    
        fstream fs( "./halo.vtk", std::fstream::out );
    
        VTK_MeshWriter_Mesh2D vtk( M, fs, "Attempt at a " );
        vtk.writeCoordinateBlock();
        vtk.writeTopDimensionalCells();

        fs.close();


    }
    
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
