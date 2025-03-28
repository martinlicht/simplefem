
#include <string>
#include <vector>

#include "../../base/include.hpp"
#include "../../utility/stl.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"



// using namespace std;

#include "vtk.testsnippet.cxx"

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: VTK output of Simplicial Mesh (2D)" << nl;
    
    {
        
        const MeshSimplicial2D Mx = StandardSquare2D_strange14();  std::string meshname = "Standard Square 2D";
//         const MeshSimplicial2D Mx = UnitedKingdom();               std::string meshname = "United Kingdom"; 
        
        internal_print( Mx, meshname );
        
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 6; c++ ) {
            
                M.uniformrefinement();

                M.shake_interior_vertices();
                
                internal_print( M, meshname, get_basename(__FILE__) );
            
            }
            
        }    
        
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 6; c++ ) {
            
                int number_of_refinement_edges = 3 + M.count_edges() / 10;
                
                std::vector<int> refinementedges;
                refinementedges.reserve( number_of_refinement_edges );
                
                for( int k = 0; k < number_of_refinement_edges; k++ )
                    refinementedges.push_back( random_integer() % M.count_edges() );
                
                sort_and_remove_duplicates( refinementedges );
                
                M.longest_edge_bisection_recursive( refinementedges );

                internal_print( M, meshname, get_basename(__FILE__) );
            
            }
            
        }    
        
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 6; c++ ) {
                
                int number_of_refinement_edges = 3 + M.count_edges() / 10;
                
                std::vector<int> refinementedges;
                refinementedges.reserve( number_of_refinement_edges );
                
                for( int k = 0; k < number_of_refinement_edges; k++ )
                    refinementedges.push_back( random_integer() % M.count_edges() );
                
                sort_and_remove_duplicates( refinementedges );
                
                M.newest_vertex_bisection_recursive( refinementedges );
                
                internal_print( M, meshname, get_basename(__FILE__) );
            
            }
        
        }
        
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
