

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;


inline void print( const MeshSimplicial2D& M, std::string meshname, const char* filename = nullptr )
{
    
    std::string basename = ( filename == nullptr ) ? getbasename(__FILE__) : filename;
    
    fstream fs( experimentfile( basename ), std::fstream::out );
    
    VTKWriter vtk( M, fs, meshname );
    vtk.writeCoordinateBlock();
    vtk.writeTopDimensionalCells();
    
    {
        FloatVector V( M.count_simplices(0), 
                        [&M](int i)->Float{
                        FloatVector point = M.getcoordinates().getvectorclone(i);
                        return std::pow( point[0], 2.0 ) + std::cos( 1.0 * 2 * 3.14159 * point[1] );
                    });
        
        vtk.writeVertexScalarData( V, "testing_scalar_data", 1.0 );
    }
    
    vtk.writeCellScalarData( FloatVector(M.count_simplices(2), 2.5 ), "cell_scalar_data" );
    
    vtk.writeCellVectorData( 
        FloatVector( M.count_simplices(2),  1.5 ),
        FloatVector( M.count_simplices(2),  2.5 ),
        FloatVector( M.count_simplices(2), -3.5 ),
        "cell_vector_data"
    );

    fs.close();
}






int main()
{
    cout << "Output of a few important meshes" << endl;
    
    
    if(true)
    {
        
        const int L = 5;
    
        MeshSimplicial2D M = SphericalSurface2D(L);
            
        cout << L << ":\t" << M.getShapemeasure() << endl;
        
        print( M, "spherical surface 2D", "sphere" );
        
    }
    
    
    
    if(true)
    {
        
        const int L = 5;
    
        MeshSimplicial2D M = LShapedDomain2D();
            
        cout << L << ":\t" << M.getShapemeasure() << endl;
        
        print( M, "L shaped domain 2D", "lshaped" );
        
    }
    
    
    
    if(true)
    {
        cout << "Unit Test for VTK output of Simplicial Mesh" << endl;
        
        const int K = 4;
        const int L = 12;
        
        MeshSimplicial2D M = Halo(K,L);
            
        print( M, "halo 2D", "halo" );

    }
    
    
    
    if(true)
    {
        const int Lmin = 3;
        const int Lmax = 9;
        
        for( int L = Lmin; L <= Lmax; L++ )
        {
            
            MeshSimplicial2D M = UnitDisk(L);
            
            cout << L << ":\t" << M.getShapemeasure() << endl;
            
            {
                fstream fs( string("./rounddisk.tex"), std::fstream::out );
                M.outputTikZ( fs );
                fs.close();
            }
            
            print( M, "round disk 2D", "rounddisk" );

        }
        
    }
    
    
    
    if(true)
    {
        
        const int Lmin = 3;
        const int Lmax = 9;
        
        MeshSimplicial2D M = Annulus( Lmin, Lmax );
        
        print( M, "annulus 2D", "annulus" );
        
    }
        
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
