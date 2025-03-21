

/**/

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main()
{
    cout << "Unit Test and Measurement for 2D NVB" << endl;
    
    {
        
        cout << "First Experiment: uniform distribution" << endl;
        
        MeshSimplicial2D M = UnitTriangle2D(); M.check();
        
        int cell_count_initial = M.count_triangles();
        int cell_marked_count  = 0;
        
        int iter_max = 4000;
        
        for( int i = 0; i < iter_max; i++ ) {
        
//             std::vector<int> markedcells;
//             markedcells.push_back( rand() % M.count_triangles() );
            
            cell_marked_count += 1;
            
            std::vector<int> markededges;
//             for( int t : markedcells ) markededges.push_back( M.get_oldest_edge( t ) );
            markededges.push_back( rand() % M.count_edges() );
            sort_and_unique( markededges );
            
//             std::cout << i << "/" << iter_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
            M.newest_vertex_bisection_recursive( markededges );
            std::cout << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
        
        }
        
        M.check();
        
    }
    
    
    {
        
        cout << "Second Experiment: fixed triangle" << endl;
        
        MeshSimplicial2D M = UnitTriangle2D(); M.check();
        
        int cell_count_initial = M.count_triangles();
        int cell_marked_count  = 0;
        
        int iter_max = 2000;
        
        for( int i = 0; i < iter_max; i++ ) {
        
//             std::vector<int> markedcells;
//             markedcells.push_back( 0 );
            
            cell_marked_count += 1;
            
            std::vector<int> markededges;
//             for( int t : markedcells ) markededges.push_back( M.get_oldest_edge( t ) );
            markededges.push_back( 0 );
            sort_and_unique( markededges );
            
            std::cout << M.get_edge_vertex( 0, 0 ) << space << M.get_edge_vertex( 0, 1 ) << nl;
            
//             std::cout << i << "/" << iter_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
            M.newest_vertex_bisection_recursive( markededges );
            std::cout << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
        
        }
        
        M.check();
        
    }
    
    
    cout << "Finished Unit Test" << endl;
    
    return 0;
}
