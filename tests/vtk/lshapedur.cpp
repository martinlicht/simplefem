

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
    
    // MeshSimplicial2D M = UnitCubeTriangulation(3,3);
    MeshSimplicial2D M = LShapedDomain2D();
    
    int l_max = 5;

    for( int l = 0; l < 5; l++ )
    {
        cout << "Print VTK-type file" << endl;
        cout << "T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << "\n";
        
        fstream fs( string("./lshaped") + std::to_string(l) + string(".vtk"), std::fstream::out );

        VTK_MeshWriter_Mesh2D vtk( M, fs, "Mein erster Test" );
        vtk.writeCoordinateBlock();
        vtk.writeTopDimensionalCells();

        fs.close();

        cout << "Refine" << endl;

        if( l != l_max ) M.uniformrefinement();
    }
        
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
