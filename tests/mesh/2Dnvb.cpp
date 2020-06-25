

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
        cout << "Unit Test for Simplicial 2D Module" << endl;
        
        {
            
            cout << "First Experiment" << endl;
            
            MeshSimplicial2D M = UnitTriangle2D();
            
            M.check();
            
//             for( int k = 0; k < 2; k++ ) M.uniformrefinement();
            
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
                sort_and_unique( markededges );
                
                std::cout << c << "/" << c_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
                M.newest_vertex_bisection_recursive( markededges );
                std::cout << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
            
            }
            
            M.check();
            
        }
        
        
        {
            
            cout << "Second Experiment" << endl;
            
            MeshSimplicial2D M = TetrahedralSurface2D(); //StandardSquare2D();
            
            M.check();
            
//             for( int k = 0; k < 2; k++ ) M.uniformrefinement();
            
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
                sort_and_unique( markededges );
                
                std::cout << c << "/" << c_max << "Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
                M.newest_vertex_bisection_recursive( markededges );
                std::cout << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
            
            }
            
            M.check();
            
        }
        
        
        {
            
            cout << "Third Experiment" << endl;
            
            MeshSimplicial2D M = TetrahedralSurface2D();
            
            M.check();
            
//             for( int k = 0; k < 2; k++ ) M.uniformrefinement();
            
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
                sort_and_unique( markededges );
                
                std::cout << c << "/" << c_max << "Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
                M.newest_vertex_bisection_recursive( markededges );
                std::cout << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
            
            }
            
            M.check();
            
        }


        cout << "Finished Unit Test" << endl;
        
        return 0;
}
