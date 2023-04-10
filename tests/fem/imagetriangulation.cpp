

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

int main()
{
    LOG << "Unit Test for Interpolation in FEEC" << nl;

    MeshSimplicial2D M = UnitSquare2D_centered();

    M.check();

    std::string name = "aurora.jpeg";
    // std::string name = "testbild.jpg";

    // PixelImage pim = readPixelImage("lena_color.tiff");
    PixelImage pim = readPixelImage( name );

    auto red   = pim.get_interpolated_red();
    auto green = pim.get_interpolated_green();
    auto blue  = pim.get_interpolated_blue();
    
    M.getcoordinates().scale( { (Float)pim.getwidth(), (Float)pim.getheight() } );

    int l_min =  0;
    int l_max =  7;
    
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

        M.uniformrefinement(); continue;

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
                samples_R[k] = red(   P[0], P[1] );
                samples_G[k] = green( P[0], P[1] );
                samples_B[k] = blue(  P[0], P[1] );
            }
            
            samples_R.shift( -interpol_red[t]   ).scale(1./K);
            samples_G.shift( -interpol_green[t] ).scale(1./K);
            samples_B.shift( -interpol_blue[t]  ).scale(1./K);

            // Float weight = flip_coin(0.9);

            Float measure = M.getMeasure( 2, t );

            Float weight =   samples_R.lpnorm(2.,measure) 
                           + samples_G.lpnorm(2.,measure) 
                           + samples_B.lpnorm(2.,measure);
            
            weights[t] = pair<int,Float>(t,weight);
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

    }

    LOG << "Finished Unit Test" << nl;

    return 0;
}
