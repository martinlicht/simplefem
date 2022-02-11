

/**/

#include <ostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../mesh/examples1D.hpp"


using namespace std;

#include "vtk.testsnippet.cxx"

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
