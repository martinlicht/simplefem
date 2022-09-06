

/**/

#include <ostream>
// #include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"

using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Simplicial 2D Module" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
        // LOG << "Unit Test for Simplicial 2D Module" << endl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        M.automatic_dirichlet_flags();

        M.check_dirichlet_flags();
            
        LOG << "Refinement" << endl;
        
        for( int c = 0; c < 5; c++ ) {

          M.midpoint_refinement_global();
          
        }
          
        
        M.check();
        
        M.check_dirichlet_flags();
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
