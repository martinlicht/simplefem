

/**/

#include "../../basic.hpp"
#include "../../utility/pixelimage.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test for Interpolation in FEEC" << nl;
        
        MeshSimplicial2D M = UnitSquare2D_strange14();

        for( int c = 0; c < 2; c++ ) M.uniformrefinement();
        
        M.check();

        PixelImage pim = readPixelImage("lena_color.tiff");

        red   = pim.get_interpolated_red();
        green = pim.get_interpolated_green();
        blue  = pim.get_interpolated_blue();
        
        FloatVector interpol_red   = Interpolation( M, M.getinnerdimension(), 2, 0, red );
        FloatVector interpol_green = Interpolation( M, M.getinnerdimension(), 2, 0, green );
        FloatVector interpol_blue  = Interpolation( M, M.getinnerdimension(), 2, 0, blue );

        fstream fs( experimentfile( getbasename(__FILE__), "svg" ), std::fstream::out );
        int num_tets = M.count_triangles();
        FloatVector red( num_tets, 128 ), green( num_tets, 240 ), blue( num_tets, 38 );
        red.random_within_range(0.,255.); green.random_within_range(0.,255.); blue.random_within_range(0.,255.);
        fs << M.outputSVG( 0.01, "array", "red", &interpol_red, &interpol_green, &interpol_blue );
        fs.close();

        
        LOG << results << nl;
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}
