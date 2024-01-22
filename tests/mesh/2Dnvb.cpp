
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
        
        for( int k = 0; k <= 2; k++ ) 
        { 
            LOG << "Uniform refinements..." << k << nl;
            M.uniformrefinement();
        }
        
        LOG << "Newest vertex bisections..." << nl;

        int cell_count_initial = M.count_triangles();
        int cell_marked_count  = 0;
        
        int c_max = 50;
        
        for( int c = 0; c < c_max; c++ ) {
        
            std::vector<int> markedcells;
            
            unsigned int p = 40;
            for( int t = 0; t < M.count_triangles(); t++ )
                if( rand() % p == 0 ) 
                    markedcells.push_back( t );
            markedcells.push_back( 0 );
            cell_marked_count += markedcells.size();
            
            std::vector<int> markededges;
            for( int t : markedcells ) markededges.push_back( M.get_oldest_edge( t ) );
            sort_and_remove_duplicates( markededges );
            
            LOG << c << "/" << c_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
            M.newest_vertex_bisection_recursive( markededges );
            LOG << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
        
        }
        
        M.check();
        
    }
    
    
    {
        
        LOG << "Second Experiment" << nl;
        
        MeshSimplicial2D M = TetrahedralSurface2D(); //StandardSquare2D();
        
        M.check();
        
        for( int k = 0; k <= 2; k++ ) 
        { 
            LOG << "Uniform refinements..." << k << nl;
            M.uniformrefinement();
        }
        
        LOG << "Newest vertex bisections..." << nl;

        int cell_count_initial = M.count_triangles();
        int cell_marked_count  = 0;
        
        int c_max = 20;
        
        for( int c = 0; c < c_max; c++ ) {
        
            std::vector<int> markedcells;
            
            unsigned int p = 10;
            for( int t = 0; t < M.count_triangles(); t++ )
                if( rand() % p == 0 ) 
                    markedcells.push_back( t );
            markedcells.push_back( 0 );
            cell_marked_count += markedcells.size();
            
            std::vector<int> markededges;
            for( int t : markedcells ) markededges.push_back( M.get_oldest_edge( t ) );
            sort_and_remove_duplicates( markededges );
            
            LOG << c << "/" << c_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
            M.newest_vertex_bisection_recursive( markededges );
            LOG << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
        
        }
        
        M.check();
        
    }
    
    
    {
        
        LOG << "Third Experiment" << nl;
        
        MeshSimplicial2D M = TetrahedralSurface2D();
        
        M.check();
        
        for( int k = 0; k <= 2; k++ ) 
        { 
            LOG << "Uniform refinements..." << k << nl;
            M.uniformrefinement();
        }
        
        LOG << "Newest vertex bisections..." << nl;

        int cell_count_initial = M.count_triangles();
        int cell_marked_count  = 0;
        
        int c_max = 50;
        
        for( int c = 0; c < c_max; c++ ) {
        
            std::vector<int> markedcells;
            
            // unsigned int p = 10;
            // for( int t = 0; t < M.count_triangles(); t++ )
            //     if( rand() % p == 0 ) 
            //         markedcells.push_back( t );
            markedcells.push_back( 0 );
            cell_marked_count += markedcells.size();
            
            std::vector<int> markededges;
            for( int t : markedcells ) markededges.push_back( M.get_oldest_edge( t ) );
            sort_and_remove_duplicates( markededges );
            
            LOG << c << "/" << c_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
            M.newest_vertex_bisection_recursive( markededges );
            LOG << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
        
        }
        
        M.check();
        
    }
        
        
    {
        
        LOG << "Fourth Experiment: uniform distribution" << nl;
        
        MeshSimplicial2D M = UnitTriangle2D(); M.check();
        
        LOG << "Newest vertex bisections..." << nl;

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
            sort_and_remove_duplicates( markededges );
            
//             LOG << i << "/" << iter_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
            M.newest_vertex_bisection_recursive( markededges );
            LOG << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
        
        }
        
        M.check();
        
    }
    
    
    {
        
        LOG << "Fifth Experiment: repeated bisection of a fixed triangle" << nl;
        
        MeshSimplicial2D M = UnitTriangle2D(); M.check();
        
        int cell_count_initial = M.count_triangles();
        int cell_marked_count  = 0;
        
        int iter_max = 2000;
        
        LOG << "Newest vertex bisections..." << nl;

        for( int i = 0; i < iter_max; i++ ) {
        
//             std::vector<int> markedcells;
//             markedcells.push_back( 0 );
            
            cell_marked_count += 1;
            
            std::vector<int> markededges;
//             for( int t : markedcells ) markededges.push_back( M.get_oldest_edge( t ) );
            markededges.push_back( 0 );
            sort_and_remove_duplicates( markededges );
            
            LOG << M.get_edge_vertex( 0, 0 ) << space << M.get_edge_vertex( 0, 1 ) << nl;
            
//             LOG << i << "/" << iter_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
            M.newest_vertex_bisection_recursive( markededges );
            LOG << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
        
        }
        
        M.check();
        
    }


        LOG << "Finished Unit Test" << nl;
        
        return 0;
}
