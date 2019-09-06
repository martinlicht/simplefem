

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
    
    const int Lmin = 3;
    const int Lmax = 6;
    
    MeshSimplicial2D M = Annulus( Lmin, Lmax );
    
    cout << "Print VTK-type file" << endl;

    fstream fs( "./annulus.vtk", std::fstream::out );

    VTK_MeshWriter_Mesh2D vtk( M, fs );
    vtk.writePreamble( "Attempt at a " );
    vtk.writeCoordinateBlock();
    vtk.writeTopDimensionalCells();
        
    fs.close();
    
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
