

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
    
    const int Lmax = 9;
    
    for( int L = 3; L <= Lmax; L++ )
    {
        
        MeshSimplicial2D M = UnitDisk(L);
        
        cout << L << ":\t" << M.getShapemeasure() << endl;
        
        {
            fstream fs( string("./rounddisk.tex"), std::fstream::out );
            M.outputTikZ( fs );
            fs.close();
        }
        
        if(false)
        {
            
            // M.uniformrefinement();
            // M.uniformrefinement();
        
            cout << "Print VTK-type file" << endl;
        
            fstream fs( string("./rounddisk") + std::to_string(L) + string(".vtk"), std::fstream::out );
        
            VTKWriter vtk( M, fs, "Attempt at a round disk" );
            vtk.writeCoordinateBlock();
            vtk.writeTopDimensionalCells();

            fs.close();

            
        }

    }
    
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
