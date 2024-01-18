
#include <cstdio>

#include <fstream>

#include "../../basic.hpp"
#include "../../utility/files.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test for Simplicial 2D Module" << nl;
        
//        MeshSimplicial2D M = TetrahedralSurface2D();
        MeshSimplicial2D M = StandardSquare2D_strange14();
        
        M.check();
        
        LOG << "Refinement" << nl;
        
        for( int c = 0; c < 2; c++ )
          M.uniformrefinement();
        
        M.check();
        
        LOG << M << nl;
        
        // std::cout << M.outputTikZ();
        puts( M.outputTikZ().c_str() );
        
        fstream fs( experimentfile( getbasename(__FILE__), "svg" ), std::fstream::out );
        int num_tets = M.count_triangles();
        FloatVector red( num_tets, 128 ), green( num_tets, 240 ), blue( num_tets, 38 );
        red.random_within_range(0.,255.); green.random_within_range(0.,255.); blue.random_within_range(0.,255.);
        fs << M.outputSVG( 0.01, "array", "white", &red, &green, &blue );
        fs.close();
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}
