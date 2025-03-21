

/**/
#include <algorithm>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "../../base/include.hpp"
#include "../../utility/files.hpp"
#include "../../utility/pixelimage.hpp"
#include "../../dense/factorization.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: Image output, San Francisco" << nl;

    /* 1. File name */

    // const std::string image_name = "blue.png";
    // const std::string image_name = "aurora.jpeg";
    // const std::string image_name = "lena_color.tiff";
    const std::string image_name = "testbild.jpg";
    // const std::string image_name = "sanfrancisco.jpg";

    const std::string execution_name = argv[0];
    const std::string name = get_parent_directory( execution_name ) + "/" + image_name;
    
    LOG << "Filename for input: " << name << nl;

    
    /* 2. Read the image and prepare the piecewise constant color functions */

    PixelImage pim = readPixelImage( name );

    auto red   = pim.get_interpolated_red();
    auto green = pim.get_interpolated_green();
    auto blue  = pim.get_interpolated_blue();
    
    
    
    MeshSimplicial2D M = UnitSquare2D_centered();

    M.check();

    // M.getCoordinates().shift( FloatVector{ 1. , 1. } );
    M.getCoordinates().scale( FloatVector{ (Float)pim.getwidth()/1., (Float)pim.getheight()/1. } );
    
    
    /* 3. Iterate over refinement levels */

    int l_min =  0;
    int l_max = 10;
    
    for( int c = 0; c < l_min; c++ ) M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ )
    {
        LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
            
        int num_volumes = M.count_triangles();

        /* 3.1 Approximate red/green/blue channels by barycentric linear functions, using linear sampling */
        
        FloatVector interpol_red(   3 * M.count_triangles() );
        FloatVector interpol_green( 3 * M.count_triangles() );
        FloatVector interpol_blue(  3 * M.count_triangles() );

        for( int t = 0; t < num_volumes; t++ ) {
                
            int K = 30;

            DenseMatrix points( K, 3 );

            FloatVector sample_R(K);
            FloatVector sample_G(K);
            FloatVector sample_B(K);

            for( int k = 0; k < K; k++ )
            {
                auto barycoords = get_random_barycentric_coordinates(2);
                assert( barycoords.is_nonnegative() );
                assert( is_numerically_one( barycoords.sum() ) );
                
                auto P = M.getPointFromBarycentric( 2, t, barycoords );

                points.setrow( k, barycoords ); // TODO(martinlicht)
                sample_R[k] = red(   P[0], P[1] );
                sample_G[k] = green( P[0], P[1] );
                sample_B[k] = blue(  P[0], P[1] );
            }
            

            auto local_red   = SolveOverconstrained( points, sample_R );
            auto local_green = SolveOverconstrained( points, sample_G );
            auto local_blue  = SolveOverconstrained( points, sample_B );
            
            // auto invpoints = Inverse(points);
            // auto local_red   = invpoints * sample_R;
            // auto local_green = invpoints * sample_G;
            // auto local_blue  = invpoints * sample_B;

            for( int c = 0; c <= 2; c++ ) {
                interpol_red[  3*t+c] = local_red  [c];
                interpol_green[3*t+c] = local_green[c];
                interpol_blue[ 3*t+c] = local_blue [c];
            }
        }


        /* 3.2 Open the output file and write color interpolation */
        
        std::fstream fs( get_available_filename( get_basename(__FILE__), "svg" ), std::fstream::out );
        fs << M.outputInterpolatingSVG( interpol_red, interpol_green, interpol_blue, 0.000, "array", "none" );
        fs.close();

        ///////////////////////////////////////////////////////////

        /* 3.3 Refinement */
        
        if( l == l_max ) break;

        // M.uniformrefinement(); continue;

        /* 3.3.1 For each triangle, we determine the weight as the integral of the interpolation from the actual values, via random sampling */
        
        std::vector<std::pair<int,Float>> weights( 
                num_volumes, std::pair<int,Float>(0,0.)
        );
        
        for( int t = 0; t < num_volumes; t++ ) {
                
            int K = 30;

            DenseMatrix points( K, 3 );

            FloatVector samples_R(K);
            FloatVector samples_G(K);
            FloatVector samples_B(K);

            for( int k = 0; k < K; k++ )
            {
                auto barycoords = get_random_barycentric_coordinates(2);
                auto P = M.getPointFromBarycentric( 2, t, barycoords );

                points.setrow( k, barycoords );
                samples_R[k] = red(   P[0], P[1] ) - barycoords[0] *   interpol_red[3*t+0] - barycoords[1] *   interpol_red[3*t+1] - barycoords[2] *   interpol_red[3*t+2];
                samples_G[k] = green( P[0], P[1] ) - barycoords[0] * interpol_green[3*t+0] - barycoords[1] * interpol_green[3*t+1] - barycoords[2] * interpol_green[3*t+2];
                samples_B[k] = blue(  P[0], P[1] ) - barycoords[0] *  interpol_blue[3*t+0] - barycoords[1] *  interpol_blue[3*t+1] - barycoords[2] *  interpol_blue[3*t+2];
            }

            Float measure = M.getMeasure( 2, t );

            Float weight =   samples_R.lpnorm(2.,measure) / K 
                           + samples_G.lpnorm(2.,measure) / K 
                           + samples_B.lpnorm(2.,measure) / K;
            // Float weight = flip_coin(0.9);
            
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
