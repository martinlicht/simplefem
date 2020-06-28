

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Simplicial 3D Module" << endl;
        
        for( int ei = 0; ei < 6; ei++ )
        {
            
            MeshSimplicial3D M = UnitSimplex3D();
            
            M.check();
            
            for( int c = 0; c < 10; c++ ) {
                M.bisect_edge( M.get_tetrahedron_edge( 0, ei ) );
            }
            
            M.check();
            
        }
        
        
        {
            
            MeshSimplicial3D M = StandardCube3D();
            
            M.check();
            
            for( int c = 0; c < 20; c++ ) {
                M.bisect_edge( M.get_tetrahedron_edge( c % 5, rand() % 6 ) );
            }
            
            M.check();
            
        }
        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
