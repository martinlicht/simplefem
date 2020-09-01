

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for VTK output of Simplicial Mesh" << endl;
    
    const int Lmin = 3;
    const int Lmax = 9;
    
    MeshSimplicial2D M = Annulus( Lmin, Lmax );
    
    cout << "Print VTK-type file" << endl;

    fstream fs( "./annulus.vtk", std::fstream::out );

    VTKWriter vtk( M, fs, "Attempt at a " );
    vtk.writeCoordinateBlock();
    vtk.writeTopDimensionalCells();
    
    {
        
        FloatVector Vx( M.count_simplices(2), 0.0 );
        FloatVector Vy( M.count_simplices(2), 1.0 );
        FloatVector Vz( M.count_simplices(2), 0.0 );
        
        vtk.writeCellVectorData( Vx, Vy, Vz, "testing_vector_data", 0.1 );
    }
                    

    
    fs.close();
    
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
