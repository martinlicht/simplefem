

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

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Simplicial 3D Module" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
        // LOG << "Unit Test for Simplicial 3D Module" << endl;
        
        {
            
            LOG << "First Experiment" << endl;
            
            MeshSimplicial3D M = UnitSimplex3D();
            
            M.check();
            
            M.automatic_dirichlet_flags();

            M.check_dirichlet_flags();

            for( int k = 0; k < 2; k++ )
                M.uniformrefinement();
            
            int cell_count_initial = M.count_tetrahedra();
            int cell_marked_count  = 0;
            
            int c_max = 3;
            
            M.check_dirichlet_flags();

            LOG << "Start iterations" << endl;
            
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
                
                LOG << c << "/" << c_max << " Refine " << markedcells.size() << "/" << M.count_tetrahedra() << " ... ";
                M.longest_edge_bisection_recursive( markededges );
                LOG << "Ratio=" << ( M.count_tetrahedra() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
                
                M.check_dirichlet_flags();

            
            
            }
            
            LOG << "check";
            
            M.check();
            
            M.check_dirichlet_flags();
            
        }
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
