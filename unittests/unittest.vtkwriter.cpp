

/**/

#include <iostream>
#include <fstream>
#include <cmath>

#include "../basic.hpp"
#include "coordinates.hpp"
#include "simplicialmesh.hpp"
#include "generatesimplicialmesh.hpp"
#include "vtkwriter.hpp"


using namespace std;

int main()
{
    
    cout << "Unit Test for VTK Writer" << endl;
    
    {
        
        SimplicialMesh M = UnitCubeTriangulation2D(5,5);
        
        cout << M << endl;

        FloatVector sampledata( M.countsimplices(0) );
        
        for( int v = 0; v < M.countsimplices(0); v++ ) {
          Float xcoord = M.getcoordinates().getdata( v, 0 );
          Float ycoord = M.getcoordinates().getdata( v, 1 );
          sampledata[v] 
          =
          std::sin( 1. * 2. * 3.14159 * xcoord) 
          *
          std::sin( 1. * 2. * 3.14159 * ycoord);
        }
        
        FloatVector samplevectorx( M.countsimplices(2) );
        FloatVector samplevectory( M.countsimplices(2) );
        FloatVector samplevectorz( M.countsimplices(2), 0. );
        
        for( int c = 0; c < M.countsimplices(2); c++ ) {
          FloatVector mid = M.getMidpoint( 2, c );
          samplevectorx[c] = mid[0];
          samplevectory[c] = mid[1];
          samplevectorz[c] = 0.; 
        }
        
        {
            
            fstream fs( "./gitter.2D.vtk", std::fstream::out );
        
            VTK_MeshWriter vtk( M, fs );
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
