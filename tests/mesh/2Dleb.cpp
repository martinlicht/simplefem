
#include <vector>

#include "../../basic.hpp"
#include "../../utility/stl.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test for Simplicial 2D Module" << nl;
        
        {
            
            LOG << "First Experiment" << nl;
            
            MeshSimplicial2D M = UnitTriangle2D();
            
            M.check();
            
            LOG << "Uniform refinements..." << nl;

            for( int k = 0; k <= 2; k++ ) {             
                LOG << "Uniform refinements..." << k << nl;
                M.uniformrefinement();
            }
            
            LOG << "Longest edge bisections..." << nl;
            
            int cell_count_initial = M.count_triangles();
            int cell_marked_count  = 0;
            
            int c_max = 6;
            
            for( int c = 0; c <= c_max; c++ ) {
            
                std::vector<int> markedcells;
                
                unsigned int p = 3;
                for( int t = 0; t < M.count_triangles(); t++ )
                    if( rand() % p == 0 ) 
                        markedcells.push_back( t );
                cell_marked_count += markedcells.size();
                
                std::vector<int> markededges;
                for( int t : markedcells ) 
                    markededges.push_back( M.get_oldest_edge( t ) );
                sort_and_remove_duplicates( markededges );
                
                LOG << c << "/" << c_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
                M.longest_edge_bisection_recursive( markededges );
                LOG << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
            
            }
            
            M.check();
            
        }
        
        
        {
            
            LOG << "Second Experiment" << nl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            LOG << "Uniform refinements..." << nl;

            for( int k = 0; k <= 2; k++ ) {             
                LOG << "Uniform refinements..." << k << nl;
                M.uniformrefinement();
            }
            
            LOG << "Longest edge bisections..." << nl;
            
            int cell_count_initial = M.count_triangles();
            int cell_marked_count  = 0;
            
            int c_max = 6;
            
            for( int c = 0; c <= c_max; c++ ) {
            
                std::vector<int> markedcells;
                
                unsigned int p = 3;
                for( int t = 0; t < M.count_triangles(); t++ )
                    if( rand() % p == 0 ) 
                        markedcells.push_back( t );
                cell_marked_count += markedcells.size();
                
                std::vector<int> markededges;
                for( int t : markedcells ) markededges.push_back( M.get_oldest_edge( t ) );
                sort_and_remove_duplicates( markededges );
                
                LOG << c << "/" << c_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
                M.longest_edge_bisection_recursive( markededges );
                LOG << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
            
            }
            
            M.check();
            
        }
        
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}
