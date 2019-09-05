

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
    
    const int Lmax = 10;
    
    for( int L = 3; L <= Lmax; L++ )
    {
        
        MeshSimplicial2D M = UnitDisk(L);
        
        cout << L << ":\t" << M.getShapemeasure() << endl;
        
        if( L == Lmax )
        {
            
            cout << "Print VTK-type file" << endl;
        
            fstream fs( "./rounddisk.vtk", std::fstream::out );
        
            VTK_MeshWriter_Mesh2D vtk( M, fs );
            vtk.writePreamble( "Attempt at a " );
            vtk.writeCoordinateBlock();
            vtk.writeTopDimensionalCells();
                
            fs.close();
            
        }

    }
    
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
