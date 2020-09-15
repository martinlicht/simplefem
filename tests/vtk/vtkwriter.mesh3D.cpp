

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
// #include "../../mesh/mesh.simplicialND.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;





inline void print( const MeshSimplicial3D& M, std::string meshname )
{
    
    fstream fs( experimentfile( getbasename(__FILE__)), std::fstream::out );
    
    VTKWriter vtk( M, fs, meshname );
    vtk.writeCoordinateBlock();
    vtk.writeTopDimensionalCells();
    
    {
        FloatVector V( M.count_simplices(0), 
                        [&M](int i)->Float{
                        FloatVector point = M.getcoordinates().getvectorclone(i);
                        return point[0] * std::sin( 10 * 2 * 3.14159 * point[1] ) + point[2];
                    });
        
        vtk.writeVertexScalarData( V, "testing_scalar_data", 1.0 );
    }
    
    vtk.writeCellScalarData( FloatVector(M.count_simplices(3), 2.5 ), "cell_scalar_data" );
    
    vtk.writeCellVectorData( 
        FloatVector( M.count_simplices(3),  1.5 ),
        FloatVector( M.count_simplices(3),  2.5 ),
        FloatVector( M.count_simplices(3), -3.5 ),
        "cell_vector_data"
    );

    fs.close();
}







int main()
{
    cout << "Unit Test for VTK output of Simplicial Mesh (3D)" << endl;
    
    {
        
        const MeshSimplicial3D Mx = StandardCube3D();  std::string meshname = "Standard Cube 3D";
        
        print( Mx, meshname );
        
        
        
        
        if(true)
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 3; c++ ) {
            
                M.uniformrefinement();
                
                print( M, meshname );
            
            }
            
        }    
        
        
        
        if(true)
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 3; c++ ) {
            
                std::vector<int> refinementedges;
                
                for( int k = 0; k < 3 + M.count_edges() / 10; k++ )
                    refinementedges.push_back( rand() % M.count_edges() );
                
                sort_and_unique( refinementedges );
                
                M.longest_edge_bisection_recursive( refinementedges );

                print( M, meshname );
            
            }
            
        }    
        

    }
    
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
