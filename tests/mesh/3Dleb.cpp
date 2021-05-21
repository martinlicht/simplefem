

/**/

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Simplicial 3D Module" << endl;
        
        {
            
            cout << "Second Experiment" << endl;
            
            MeshSimplicial3D M = UnitSimplex3D();
            
            M.check();
            
            for( int k = 0; k < 2; k++ ) M.uniformrefinement();
            
            int cell_count_initial = M.count_tetrahedra();
            int cell_marked_count  = 0;
            
            int c_max = 6;
            
            for( int c = 0; c < c_max; c++ ) {
            
                std::vector<int> markedcells;
                
                unsigned int p = 3;
                for( int t = 0; t < M.count_tetrahedra(); t++ )
                    if( rand() % p == 0 ) 
                        markedcells.push_back( t );
                cell_marked_count += markedcells.size();
                
                std::vector<int> markededges;
                for( int t : markedcells ) markededges.push_back( M.get_oldest_edge( t ) );
                sort_and_remove_duplicates( markededges );
                
                std::cout << c << "/" << c_max << " Refine " << markedcells.size() << "/" << M.count_tetrahedra() << " ... ";
                M.longest_edge_bisection_recursive( markededges );
                std::cout << "Ratio=" << ( M.count_tetrahedra() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
            
            }
            
            M.check();
            
        }
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
