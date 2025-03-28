

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

int main()
{
    clog << "Unit Test for VTK output of Simplicial Mesh" << endl;
    
    {
        
        // MeshSimplicial3D M = UnitCubeTriangulation(3,3);
        MeshSimplicial3D M = UnitCube3D();

        clog << "Conduct refinements" << nl;
        
        for( int c = 0; c < 2; c++ ) {
        
            M.uniformrefinement();
        
        }
        

        

        
        clog << "Print VTK-type file" << endl;
        
        fstream fs( "./gitter.vtk", std::fstream::out );
        
        {
            
            VTKWriter vtk( M, fs, "Mein erster Test" );
            vtk.writeCoordinateBlock();
            vtk.writeTopDimensionalCells();
            
        }
        
        fs.close();

        }
    
        
    
    clog << "Finished Unit Test" << endl;

    return 0;
}
