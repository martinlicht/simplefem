

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../vtk/vtkwriter.hpp"
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
            
            VTKWriter vtk( M, fs, "One-dimensional Test Mesh" );
            vtk.writeCoordinateBlock();
            vtk.writeTopDimensionalCells();
            
            vtk.writeVertexScalarData( FloatVector(M.count_simplices(0), 1.5 ), "vertex_scalar_data" );
            vtk.writeCellScalarData( FloatVector(M.count_simplices(1), 2.5 ), "cell_scalar_data" );
            vtk.writeCellVectorData( 
                FloatVector(M.count_simplices(1),  1.5 ),
                FloatVector(M.count_simplices(1),  2.5 ),
                FloatVector(M.count_simplices(1), -3.5 ),
                "cell_vector_data"
            );
            
            
            
        }
        
        fs.close();
        
    }
    
    cout << "Finished Unit Test" << endl;
    
    return 0;
}
