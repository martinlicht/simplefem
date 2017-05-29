

/**/

#include <iostream>
#include <fstream>
#include <cmath>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/manifold.2D.hpp"
#include "../mesh/vtkwriter.manifold2D.hpp"


using namespace std;

int main()
{
    
    cout << "Unit Test for VTK Writer on 2D manifolds" << endl;
    
    {
        
        ManifoldTriangulation2D M = UnitSquare();
        
        cout << M << endl;

        FloatVector sampledata( M.count_vertices() );
        
        assert( M.count_vertices() == M.getcoordinates().getnumber() );
        
        for( int v = 0; v < M.count_vertices(); v++ ) {
          Float xcoord = M.getcoordinates().getdata( v, 0 );
          Float ycoord = M.getcoordinates().getdata( v, 1 );
          sampledata[v] 
          =
          std::sin( 1. * 2. * 3.14159 * xcoord) 
          *
          std::sin( 1. * 2. * 3.14159 * ycoord) + 0.1;
        }
        
        FloatVector samplevectorx( M.count_triangles() );
        FloatVector samplevectory( M.count_triangles() );
        FloatVector samplevectorz( M.count_triangles(), 0. );
        
        for( int t = 0; t < M.count_triangles(); t++ ) {
          FloatVector mid = M.get_triangle_midpoint( t );
//           Float mid[2] = { 0.1234, 5.6789 };
          samplevectorx[t] = mid[0];
          samplevectory[t] = mid[1];
          samplevectorz[t] = 0.; 
        }
        
        {
            
            fstream fs( "./gitter.2D.manifold.vtk", std::fstream::out );
        
            VTK_MeshWriter_Manifold2D vtk( M, fs );
            vtk.writePreamble( "Two-dimensional Test" );
            vtk.writeCoordinateBlock();
            vtk.writeTopDimensionalCells();
        
            vtk.writeVertexScalarData( sampledata, "Beispieldaten" );
        
            vtk.writeCellVectorData( samplevectorx, samplevectory, samplevectorz, "Beispielvektor" );
        
            fs.close();

        }
        
        
    }
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
