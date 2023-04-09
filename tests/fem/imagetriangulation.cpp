

/**/

#include "../../basic.hpp"
#include "../../utility/files.hpp"
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

        for( int c = 0; c < 6; c++ ) M.uniformrefinement();
        
        M.check();

        PixelImage pim = readPixelImage("lena_color.tiff");

        auto red   = pim.get_interpolated_red();
        auto green = pim.get_interpolated_green();
        auto blue  = pim.get_interpolated_blue();
        
        FloatVector interpol_red   = Interpolation( M, M.getinnerdimension(), 0, 0, [&]( const FloatVector& vec ){ return FloatVector({ red  (vec[0],vec[1]) }); } );
        FloatVector interpol_green = Interpolation( M, M.getinnerdimension(), 0, 0, [&]( const FloatVector& vec ){ return FloatVector({ green(vec[0],vec[1]) }); } );
        FloatVector interpol_blue  = Interpolation( M, M.getinnerdimension(), 0, 0, [&]( const FloatVector& vec ){ return FloatVector({ blue (vec[0],vec[1]) }); } );

        fstream fs( experimentfile( getbasename(__FILE__), "svg" ), std::fstream::out );
        int num_tets = M.count_triangles();
        
        interpol_red.to_absolute();
        interpol_green.to_absolute();
        interpol_blue.to_absolute();
        
        // FloatVector redvec( num_tets, 128 ), greenvec( num_tets, 240 ), bluevec( num_tets, 38 );
        // redvec.random_within_range(0.,255.); greenvec.random_within_range(0.,255.); blue.random_within_range(0.,255.);
        fs << M.outputSVG( 0.000, "array", "red", &interpol_red, &interpol_green, &interpol_blue );
        fs.close();

        
        LOG << interpol_red.min() << space << interpol_red.maxnorm() << nl;

        LOG << "Finished Unit Test" << nl;
        
        return 0;
}
