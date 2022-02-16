

/**/

#include <ostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;


#include "vtk.testsnippet.cxx"






int main()
{
    LOG << "Output of a few important meshes" << endl;
    
    
    //if(false)
    {
        
        const int L = 5;
    
        MeshSimplicial2D M = SphericalSurface2D(L);
            
        LOG << L << ":\t" << M.getShapemeasure() << endl;
        
        internal_print( M, "spherical surface 2D", "sphere" );
        
    }
    
    
    
    //if(false)
    {
        
        const int L = 5;
    
        MeshSimplicial2D M = LShapedDomain2D();
            
        LOG << L << ":\t" << M.getShapemeasure() << endl;
        
        internal_print( M, "L shaped domain 2D", "lshaped" );
        
    }
    
    
    
    //if(false)
    {
        LOG << "Unit Test for VTK output of Simplicial Mesh" << endl;
        
        const int K = 4;
        const int L = 12;
        
        MeshSimplicial2D M = Halo(K,L);
            
        internal_print( M, "halo 2D", "halo" );

    }
    
    
    
    //if(false)
    {
        const int Lmin = 3;
        const int Lmax = 9;
        
        for( int L = Lmin; L <= Lmax; L++ )
        {
            
            MeshSimplicial2D M = UnitDisk(L);
            
            LOG << L << ":\t" << M.getShapemeasure() << endl;
            
            {
                fstream fs( string("./rounddisk.tex"), std::fstream::out );
                M.outputTikZ( fs );
                fs.close();
            }
            
            internal_print( M, "round disk 2D", "rounddisk" );

        }
        
    }
    
    
    
    //if(false)
    {
        
        const int Lmin = 3;
        const int Lmax = 9;
        
        MeshSimplicial2D M = Annulus( Lmin, Lmax );
        
        internal_print( M, "annulus 2D", "annulus" );
        
    }
        
    
    LOG << "Finished Unit Test" << endl;

    return 0;
}
