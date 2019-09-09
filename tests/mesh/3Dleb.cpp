

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
            
            MeshSimplicial3D M = UnitSimplex3D();
            
            M.check();
            
            for( int c = 0; c < 1; c++ ) {
            
                std::vector<int> refinementedges;
                
                for( int k = 0; k < 3 + M.count_edges() / 10; k++ )
                    refinementedges.push_back( rand() % M.count_edges() );
                
                std::sort( refinementedges.begin(), refinementedges.end() );
                auto last = std::unique( refinementedges.begin(), refinementedges.end() );
                refinementedges.erase( last, refinementedges.end() );
                
                M.longest_edge_bisection( refinementedges );
            
            }
            
            M.check();
            
        }
        
        
        {
            
            MeshSimplicial3D M = StandardCube3D();
            
            M.check();
            
            for( int c = 0; c < 1; c++ ) {
            
                std::vector<int> refinementedges;
                
                for( int k = 0; k < 3 + M.count_edges() / 10; k++ )
                    refinementedges.push_back( rand() % M.count_edges() );
                
                std::sort( refinementedges.begin(), refinementedges.end() );
                auto last = std::unique( refinementedges.begin(), refinementedges.end() );
                refinementedges.erase( last, refinementedges.end() );
                
                M.longest_edge_bisection( refinementedges );
            
            }
            
            M.check();
            
        }
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
