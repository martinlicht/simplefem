

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
    
    const int L = 5;
    
    MeshSimplicial2D M = SphericalSurface2D(L);
        
    cout << L << ":\t" << M.getShapemeasure() << endl;
        
            
            cout << "Print VTK-type file" << endl;
        
            fstream fs( "./sphere.vtk", std::fstream::out );
        
            VTKWriter vtk( M, fs, "Attempt at a spherical surface" );
            vtk.writeCoordinateBlock();
            vtk.writeTopDimensionalCells();
                
            {
                FloatVector V( M.count_simplices(2), 
                                [&M](int i)->Float{
                                    FloatVector point = M.get_triangle_midpoint(i);
                                    return std::pow( point[0], 2.0 ) + std::cos( 1.0 * 2 * 3.14159 * point[1] );
                            });
                
                vtk.writeCellScalarData(V,"testing_scalar_data",1.0);
            }

            fs.close();

            
    
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
