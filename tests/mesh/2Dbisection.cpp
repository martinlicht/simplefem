
#include "../../basic.hpp"
#include "../../utility/random.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Simplicial 2D Module" << nl;
    
    for( int ei = 0; ei < 3; ei++ )
    {
        
        MeshSimplicial2D M = UnitTriangle2D();
        
        M.check();
        
        for( int c = 0; c < 10; c++ ) {
            M.bisect_edge( M.get_triangle_edge( 0, ei ) );
        }
        
        M.check();
        
    }
    
    
    {
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        for( int c = 0; c < 20; c++ ) {
            M.bisect_edge( M.get_triangle_edge( c % 2, random_integer() % 3 ) );
        }
        
        M.check();
        
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
