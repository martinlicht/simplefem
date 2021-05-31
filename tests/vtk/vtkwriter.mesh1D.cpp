

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../operators/floatvector.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../mesh/examples1D.hpp"


using namespace std;






inline void internal_print( const MeshSimplicial1D& M, std::string meshname )
{
    
    fstream fs( experimentfile( getbasename(__FILE__)), std::fstream::out );
    
    VTKWriter vtk( M, fs, meshname );
    vtk.writeCoordinateBlock();
    vtk.writeTopDimensionalCells();
    
    {
        FloatVector V( M.count_simplices(0), 
                        [&M](int i)->Float{
                        FloatVector point = M.getcoordinates().getvectorclone(i);
                        return std::sin( 10 * 2 * 3.14159 * point[1] );
                    });
        
        vtk.writeVertexScalarData( V, "testing_scalar_data", 1.0 );
    }
    
    vtk.writeCellScalarData( FloatVector(M.count_simplices(1), 2.5 ), "cell_scalar_data" );
    
    vtk.writeCellVectorData( 
        FloatVector( M.count_simplices(1),  1.5 ),
        FloatVector( M.count_simplices(1),  2.5 ),
        FloatVector( M.count_simplices(1), -3.5 ),
        "cell_vector_data"
    );

    fs.close();
}







int main()
{
    LOG << "Unit Test for VTK output of Simplicial Mesh (1D)" << endl;
    
    {
        
        MeshSimplicial1D Mx = StandardInterval1D(); string meshname = string("One-dimensional Test Mesh: ") + getbasename(__FILE__);
        
        internal_print( Mx, meshname );
        
        
        
        
        //if(false)
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 6; c++ ) {
            
                M.uniformrefinement();
                
                internal_print( M, meshname );
            
            }
            
        }    
        
        
        
        
        
        
        

    }
    
    LOG << "Finished Unit Test" << endl;
    
    return 0;
}
