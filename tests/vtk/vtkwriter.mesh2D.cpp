

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
    
    {
        
        // MeshSimplicial2D M = UnitCubeTriangulation(3,3);
        MeshSimplicial2D M = UnitSquare2D();
        
        cout << M << endl;

                cout << "Print VTK-type file" << endl;
                
                fstream fs( "./gitter.vtk", std::fstream::out );
                
                {
                    
                    VTKWriter vtk( M, fs, "Mein erster Test" );
                    vtk.writeCoordinateBlock();
                    vtk.writeTopDimensionalCells();
                    
                    {
                        FloatVector V( M.count_simplices(0), 
                                     [&M](int i)->Float{
                                        FloatVector point = M.getcoordinates().getvectorclone(i);
                                        return point[0] * std::sin( 10 * 2 * 3.14159 * point[1] );
                                    });
                        
                        vtk.writeVertexScalarData(V,"testing_scalar_data",1.0);
                    }
                    
                }
                
                fs.close();

        }
    
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
