
#include <vector>

#include "../../basic.hpp"
#include "../../utility/stl.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test for Simplicial 3D Module" << nl;
        
        {
            
//             MeshSimplicial3D M = UnitSimplex3D();
            MeshSimplicial3D M = StandardCubeFive3D();
            
            M.check();
            
            LOG << "Uniform refinements..." << nl;

            for( int k = 0; k <= 2; k++ ) {             
                LOG << "Uniform refinements..." << k << nl;
                M.uniformrefinement();
            }
            
            LOG << "Longest edge bisections..." << nl;

            int cell_count_initial = M.count_tetrahedra();
            int cell_marked_count  = 0;
            
            int c_max = 4;
            
            for( int c = 0; c <= c_max; c++ ) {
            
                std::vector<int> markedcells;
                
                unsigned int p = 3;
                for( int t = 0; t < M.count_tetrahedra(); t++ )
                    if( rand() % p == 0 ) 
                        markedcells.push_back( t );
                cell_marked_count += markedcells.size();
                
                std::vector<int> markededges;
                for( int t : markedcells ) markededges.push_back( M.get_oldest_edge( t ) );
                sort_and_remove_duplicates( markededges );
                
                LOG << c << "/" << c_max << " Refine " << markedcells.size() << "/" << M.count_tetrahedra() << " ... ";
                M.longest_edge_bisection_recursive( markededges );
                LOG << "Ratio=" << ( M.count_tetrahedra() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
            
            }
            
            M.check();
            
        }
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}
