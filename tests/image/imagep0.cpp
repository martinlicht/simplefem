

/**/
#include <algorithm>
#include <vector>

#include "../../basic.hpp"
#include "../../utility/files.hpp"
#include "../../utility/pixelimage.hpp"
#include "../../utility/random.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Image triangulation" << nl;
    
    MeshSimplicial2D M = UnitSquare2D_centered();

    M.check();

//     const std::string image_name = "aurora.jpeg";
//     const std::string image_name = "lena_color.tiff";
//     const std::string image_name = "testbild.jpg";
    const std::string image_name = "sanfrancisco.jpg";
    
    const std::string execution_name = argv[0];
    const std::string name = get_parent_directory( execution_name ) + "/" + image_name;
    
    LOG << name << nl;

    // PixelImage pim = readPixelImage("lena_color.tiff");
    PixelImage pim = readPixelImage( name );

    auto red   = pim.get_interpolated_red();
    auto green = pim.get_interpolated_green();
    auto blue  = pim.get_interpolated_blue();
    
    M.getcoordinates().scale( 
        FloatVector({ 
            (Float)pim.getwidth(), (Float)pim.getheight() 
        }) 
    );

    int l_min =  0;
    int l_max =  10;
    
    for( int c = 0; c < l_min; c++ ) M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ )
    {
        LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
            
        int num_volumes = M.count_triangles();

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
                // auto P = M.get_midpoint(2,t);
                sample_R += red(   P[0], P[1] ) / K;
                sample_G += green( P[0], P[1] ) / K;
                sample_B += blue(  P[0], P[1] ) / K;
            }
            
            interpol_red[t]   = sample_R;
            interpol_green[t] = sample_G;
            interpol_blue[t]  = sample_B;
        }


        // FloatVector redvec( num_tets, 128 ), greenvec( num_tets, 240 ), bluevec( num_tets, 38 );
        // redvec.random_within_range(0.,255.); greenvec.random_within_range(0.,255.); blue.random_within_range(0.,255.);

        fstream fs( experimentfile( getbasename(__FILE__), "svg" ), std::fstream::out );
        fs << M.outputSVG( 0.000, "array", "none", &interpol_red, &interpol_green, &interpol_blue );
        fs.close();

        if( l == l_max ) break;

        // M.uniformrefinement(); continue;

        ///////////////////////////////////////////////////////////

        // container for all the weights 
        vector<pair<int,Float>> weights( 
                num_volumes, std::pair<int,Float>(0,0.)
        );
        
        // compute the weight of all the volumes
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
            // Float weight = flip_coin(0.9);
            
            weights[t] = pair<int,Float>( t, weight );
        }

        // sort descending by weight
        std::sort( weights.begin(), weights.end(), 
            []( const pair<int,Float>& a, const pair<int,Float>& b )
            { return a.second > b.second; }
        );
        
        // list all indices to be refined  
        std::vector<int> to_refine;
        int share = 3;
        LOG << weights.front().second << space << weights.back().second << nl;
        to_refine.reserve( num_volumes / share + 1 );
        for( int i = 0; i < num_volumes / share + 1; i++ )
            to_refine.push_back( weights[i].first );

        // dew it!
        M.longest_edge_bisection_recursive( to_refine );

        ///////////////////////////////////////////////////////////

    }

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
