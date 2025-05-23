

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
    clog << "Unit Test for VTK output of Simplicial Mesh" << endl;
    
    {
        
        // MeshSimplicial2D M = UnitCubeTriangulation(3,3);
        MeshSimplicial2D M = LShapedDomain2D();

        clog << "Conduct refinements" << nl;
        
        for( int c = 0; c < 4; c++ ) {
        
            M.uniformrefinement();
        
        }
        
        Coordinates& C = M.getcoordinates();
        for( int c = 0; c < C.getnumber(); c++ )
            C.setdata( c, 2, 2. * C.getdata(c,0) + 3. * C.getdata(c,1) );

        
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
