

/**/

#include <ostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"



using namespace std;

#include "vtk.testsnippet.cxx"

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "VTK output of Simplicial Mesh (2D)" );

int main()
{
    LOG << "Unit Test: " << TestName << endl;
    
    {
        
        const MeshSimplicial2D Mx = StandardSquare2D_strange14();  std::string meshname = "Standard Square 2D";
//         const MeshSimplicial2D Mx = UnitedKingdom();               std::string meshname = "United Kingdom"; 
        
        internal_print( Mx, meshname );
        
        
        
        
        //if(false)
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 6; c++ ) {
            
                M.uniformrefinement();
                
                internal_print( M, meshname );
            
            }
            
        }    
        
        
        
        //if(false)
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 6; c++ ) {
            
                std::vector<int> refinementedges;
                
                for( int k = 0; k < 3 + M.count_edges() / 10; k++ )
                    refinementedges.push_back( rand() % M.count_edges() );
                
                sort_and_remove_duplicates( refinementedges );
                
                M.longest_edge_bisection_recursive( refinementedges );

                internal_print( M, meshname );
            
            }
            
        }    
        
        
        
        //if(false)
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 6; c++ ) {
                
                std::vector<int> refinementedges;
                
                for( int k = 0; k < 3 + M.count_edges() / 10; k++ )
                    refinementedges.push_back( rand() % M.count_edges() );
                
                sort_and_remove_duplicates( refinementedges );
                
                M.newest_vertex_bisection_recursive( refinementedges );
                
                internal_print( M, meshname );
            
            }
        
        }
        
        
        

    }
    
        
    
    LOG << "Finished Unit Test: " << TestName << endl;

    return 0;
}
