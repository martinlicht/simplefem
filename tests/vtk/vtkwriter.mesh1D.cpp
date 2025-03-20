
#include <string>

#include "../../base/include.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"


// using namespace std;

#include "vtk.testsnippet.cxx"

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: VTK output of Simplicial Mesh (1D)" << nl;
    
    {
        
        MeshSimplicial1D Mx = StandardInterval1D(); std::string meshname = std::string("One-dimensional Test Mesh: ") + get_basename(__FILE__);
        
        internal_print( Mx, meshname );
        
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 6; c++ ) {
            
                M.uniformrefinement();

                // M.shake_interior_vertices(); // only if inner and outer dimension match
                
                internal_print( M, meshname );
            
            }
            
        }    
        
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
