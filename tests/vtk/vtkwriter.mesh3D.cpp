

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.hpp"


using namespace std;

#include "vtk.testsnippet.cxx"

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "VTK output of Simplicial Mesh (3D)" );

int main()
{
    LOG << "Unit Test: " << TestName << endl;
    
    {
        
        const MeshSimplicial3D Mx = StandardCube3D();  std::string meshname = "Standard Cube 3D";
        
        internal_print( Mx, meshname );
        
        
        
        
        //if(false)
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 3; c++ ) {
            
                M.uniformrefinement();
                
                internal_print( M, meshname );
            
            }
            
        }    
        
        
        
        //if(false)
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 3; c++ ) {
            
                std::vector<int> refinementedges;
                
                for( int k = 0; k < 3 + M.count_edges() / 10; k++ )
                    refinementedges.push_back( rand() % M.count_edges() );
                
                sort_and_remove_duplicates( refinementedges );
                
                M.longest_edge_bisection_recursive( refinementedges );

                internal_print( M, meshname );
            
            }
            
        }    
        

    }
    
        
    
    LOG << "Finished Unit Test: " << TestName << endl;

    return 0;
}
