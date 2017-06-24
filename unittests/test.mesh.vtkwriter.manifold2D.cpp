

/**/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/mesh.manifold2D.hpp"
#include "../mesh/vtkwriter.manifold2D.hpp"


using namespace std;

int main()
{
    
    cout << "Unit Test for VTK Writer on 2D manifolds" << endl;
    
    {
        
        MeshManifold2D M =  UnitSquare();
        
        cout << M << endl;

        srand( 4711 );
        
//         for( int c = 0; c < 10; c++ ) 
//         for( int c = 0; c < 100; c++ ) 
//         {
//           int e = 7 + c % 2; c % M.count_edges();
//           cout << "Bisect edge: " << e << nl;
//           
//           M.bisect_edge( e );
//           
// //           cout << M << endl;
//           
//           M.check();
//         }
        
        for( int c = 0; c < 6; c++ )
        {
          cout << "Uniform refinement: " << c << endl;
          M.uniformrefinement();
          
          cout << M << endl;
          cout << "check " << c << endl;
          M.check();
          
        }
        
        FloatVector sampledata( M.count_vertices() );
        
        assert( M.count_vertices() == M.getcoordinates().getnumber() );
        
        for( int v = 0; v < M.count_vertices(); v++ ) {
          Float xcoord = M.getcoordinates().getdata( v, 0 );
          Float ycoord = M.getcoordinates().getdata( v, 1 );
          sampledata[v] 
          =
          std::sin( 1. * 2. * 3.14159 * xcoord) 
          *
          std::sin( 1. * 2. * 3.14159 * ycoord) + 0.0;
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
