

/**/
#include <algorithm>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "../../base/include.hpp"
#include "../../utility/files.hpp"
#include "../../utility/pixelimage.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: Image triangulation" << nl;
    
    MeshSimplicial2D M = UnitSquare2D_centered();

    M.check();

    /* 1. File name */

    // const std::string image_name = "aurora.jpeg";
    // const std::string image_name = "lena_color.tiff";
    // const std::string image_name = "testbild.jpg";
    const std::string image_name = "sanfrancisco.jpg";
    
    const std::string execution_name = argv[0];
    const std::string name = get_parent_directory( execution_name ) + "/" + image_name;
    
    LOG << "Filename for input: " << name << nl;

    
    /* 2. Read the image and prepare the piecewise constant color functions */

    PixelImage pim = readPixelImage( name );

    auto red   = pim.get_interpolated_red();
    auto green = pim.get_interpolated_green();
    auto blue  = pim.get_interpolated_blue();
    
    M.getCoordinates().scale( 
        FloatVector({ 
            (Float)pim.getwidth(), (Float)pim.getheight() 
        }) 
    );

    
    /* 3. Iterate over refinement levels */

    int l_min =  0;
    int l_max =  10;
    
    for( int c = 0; c < l_min; c++ ) M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ )
    {
        LOG << "Level: " << l_min << " <= " << l << " <= " << l_max << nl;
            
        int num_volumes = M.count_triangles();

        /* 3.1 Approximate the average red/green/blue value of each triangle via random sampling */
        
        FloatVector interpol_red(   M.count_triangles() );
        FloatVector interpol_green( M.count_triangles() );
        FloatVector interpol_blue(  M.count_triangles() );

        // interpol_red   = Interpolation( M, M.getinnerdimension(), 0, 0, [&]( const FloatVector& vec ){ return FloatVector({ red  (vec[0],vec[1]) }); } );
        // interpol_green = Interpolation( M, M.getinnerdimension(), 0, 0, [&]( const FloatVector& vec ){ return FloatVector({ green(vec[0],vec[1]) }); } );
        // interpol_blue  = Interpolation( M, M.getinnerdimension(), 0, 0, [&]( const FloatVector& vec ){ return FloatVector({ blue (vec[0],vec[1]) }); } );

        for( int t = 0; t < num_volumes; t++ ) {
                
            int K = 30;
            Float sample_R = 0.;
            Float sample_G = 0.;
            Float sample_B = 0.;

            for( int k = 0; k < K; k++ )
            {
                auto P = M.get_random_point(2,t);
                sample_R += red(   P[0], P[1] );
                sample_G += green( P[0], P[1] );
                sample_B += blue(  P[0], P[1] );
            }
            
            interpol_red[t]   = sample_R / K;
            interpol_green[t] = sample_G / K;
            interpol_blue[t]  = sample_B / K;
        }

        /* 3.2 Open the output file and write color interpolation */
        
        std::fstream fs( get_available_filename( get_basename(__FILE__), "svg" ), std::fstream::out );
        fs << M.outputSVG( 0.000, "array", "none", &interpol_red, &interpol_green, &interpol_blue );
        fs.close();

        ///////////////////////////////////////////////////////////

        /* 3.3 Refinement */
        
        if( l == l_max ) break;

        // M.uniformrefinement(); continue;

        /* 3.3.1 For each triangle, we determine the weight as the integral of the deviation from the mean */
        
        std::vector<std::pair<int,Float>> weights( 
                num_volumes, std::pair<int,Float>(0,0.)
        );
        
        for( int t = 0; t < num_volumes; t++ ) {
                
            int K = 30;
            FloatVector samples_R(K,0.);
            FloatVector samples_G(K,0.);
            FloatVector samples_B(K,0.);

            for( int k = 0; k < K; k++ )
            {
                auto P = M.get_random_point(2,t);
                samples_R[k] = red(   P[0], P[1] ) - interpol_red[t];
                samples_G[k] = green( P[0], P[1] ) - interpol_green[t];
                samples_B[k] = blue(  P[0], P[1] ) - interpol_blue[t];
            }
            
            Float measure = M.getMeasure( 2, t );

            Float weight =   samples_R.lpnorm(2.,measure) / K 
                           + samples_G.lpnorm(2.,measure) / K 
                           + samples_B.lpnorm(2.,measure) / K;
            
            weights[t] = std::pair<int,Float>( t, weight );
        }

        /* 3.3.2 Sort the triangles by their weight */
        
        std::sort( weights.begin(), weights.end(), 
            []( const std::pair<int,Float>& a, const std::pair<int,Float>& b )
            { return a.second > b.second; }
        );
        
        LOG << "Weights min/max: " << weights.front().second << space << weights.back().second << nl;
        
        /* 3.3.2 Refine the upper third of the ranked triangles */
        
        int share = 3;
        std::vector<int> to_refine;
        to_refine.reserve( num_volumes / share + 1 );
        
        for( int i = 0; i < num_volumes / share + 1; i++ )
            to_refine.push_back( weights[i].first );

        M.longest_edge_bisection_recursive( to_refine );

        ///////////////////////////////////////////////////////////

    }

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
