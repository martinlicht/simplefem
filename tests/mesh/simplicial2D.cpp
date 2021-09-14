

/**/

#include <iostream>
#include <fstream>

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
        
//        MeshSimplicial2D M = TetrahedralSurface2D();
        MeshSimplicial2D M = StandardSquare2D_strange14();
        
        M.check();
        
        LOG << "Refinement" << endl;
        
        for( int c = 0; c < 2; c++ )
          M.uniformrefinement();
        
        M.check();
        
        LOG << M << endl;
        
        M.outputTikZ( std::cout );
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
