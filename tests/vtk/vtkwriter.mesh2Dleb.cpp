

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../vtk/vtkwriter.mesh2D.hpp"


using namespace std;

int main()
{
    clog << "Unit Test for VTK output of Simplicial Mesh" << endl;
    
    {
        
        // MeshSimplicial2D M = UnitCubeTriangulation(3,3);
        MeshSimplicial2D M = LShapedDomain2D();

        clog << "Conduct refinements" << nl;
        
        for( int c = 0; c < 4; c++ ) {
        
            std::vector<int> refinementedges;
            
            for( int k = 0; k < 3 + M.count_edges() / 10; k++ )
                refinementedges.push_back( rand() % M.count_edges() );
            
            std::sort( refinementedges.begin(), refinementedges.end() );
            auto last = std::unique( refinementedges.begin(), refinementedges.end() );
            refinementedges.erase( last, refinementedges.end() );
            
            M.longest_edge_bisection( refinementedges );
        
        }
        
        
        
        Coordinates& C = M.getcoordinates();
        for( int c = 0; c < C.getnumber(); c++ )
            C.setdata( c, 2, 2. * C.getdata(c,0) + 3. * C.getdata(c,1) );

        
        
        clog << "Print VTK-type file" << endl;
        
        fstream fs( "./gitter.vtk", std::fstream::out );
        
        {
            
            VTK_MeshWriter_Mesh2D vtk( M, fs );
            vtk.writePreamble( "Mein erster Test" );
            vtk.writeCoordinateBlock();
            vtk.writeTopDimensionalCells();
            
        }
        
        fs.close();

        }
    
        
    
    clog << "Finished Unit Test" << endl;

    return 0;
}
