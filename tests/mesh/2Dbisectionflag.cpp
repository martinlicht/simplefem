

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test for Simplicial 2D Module";// << endl;
        
        for( int ei = 0; ei < 3; ei++ )
        {
            
            MeshSimplicial2D M = UnitTriangle2D();
            
            M.check();

            M.automatic_dirichlet_flags();
            
            M.check_dirichlet_flags();
            
            for( int c = 0; c < 10; c++ ) {
                M.bisect_edge( M.get_triangle_edge( 0, ei ) );
            }
            
            M.check();
            
            M.check_dirichlet_flags();
            
        }
        
        
        {
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();

            M.automatic_dirichlet_flags();
            
            M.check_dirichlet_flags();
            
            for( int c = 0; c < 10; c++ ) {
                M.bisect_edge( M.get_triangle_edge( c % 2, rand() % 3 ) );
            }
            
            M.check();
            
            M.check_dirichlet_flags();
            
        }
        
        LOG << "Finished Unit Test";// << endl;
        
        return 0;
}
