

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
    
    const int Lmax = 6;
    
    for( int L = 3; L <= Lmax; L++ )
    {
        
        MeshSimplicial2D M = UnitDisk(L);
        
        cout << L << ":\t" << M.getShapemeasure() << endl;
        
        if( L == Lmax )
        {
            
            M.uniformrefinement();
            M.uniformrefinement();
        
            cout << "Print VTK-type file" << endl;
        
            fstream fs( "./rounddisk.vtk", std::fstream::out );
        
            VTK_MeshWriter_Mesh2D vtk( M, fs );
            vtk.writePreamble( "Attempt at a " );
            vtk.writeCoordinateBlock();
            vtk.writeTopDimensionalCells();
                
            {
                FloatVector V( M.count_simplices(0), 
                                [&M](int i)->Float{
                                FloatVector point = M.getcoordinates().getvectorclone(i);
                                return point[0] + std::cos( 2.0 * 2 * 3.14159 * point[1] );
                            });
                
                vtk.writeVertexScalarData(V,"testing_scalar_data",1.0);
            }

            fs.close();

            
        }

    }
    
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
