

/**/
#include <algorithm>
#include <vector>

#include "../../basic.hpp"
#include "../../utility/files.hpp"
#include "../../utility/pixelimage.hpp"
#include "../../utility/random.hpp"
#include "../../dense/functions.hpp"
#include "../../dense/factorization.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Interpolation in FEEC" << nl;

    MeshSimplicial2D M = UnitSquare2D_centered();

    M.check();

    // const std::string image_name = "blue.png";
    // const std::string image_name = "aurora.jpeg";
    // const std::string image_name = "lena_color.tiff";
    // const std::string image_name = "testbild.jpg";
    std::string image_name = "sanfrancisco.jpg";

    const std::string execution_name = argv[0];
    const std::string name = get_parent_directory( execution_name ) + "/" + image_name;
    
    LOG << name << nl;
    
    PixelImage pim = readPixelImage( name );

    auto red   = pim.get_interpolated_red();
    auto green = pim.get_interpolated_green();
    auto blue  = pim.get_interpolated_blue();
    
    M.getcoordinates().scale( 
        FloatVector{ 
            (Float)pim.getwidth(), (Float)pim.getheight() 
        }
    );
    
    int l_min =  0;
    int l_max = 10;
    
    for( int c = 0; c < l_min; c++ ) M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ )
    {
        LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
            
        int num_volumes = M.count_triangles();

        FloatVector interpol_red(   3 * M.count_triangles() );
        FloatVector interpol_green( 3 * M.count_triangles() );
        FloatVector interpol_blue(  3 * M.count_triangles() );

        for( int t = 0; t < num_volumes; t++ ) {
                
            int K = 30;

            DenseMatrix points( K, 3 );

            FloatVector sample_R(K);
            FloatVector sample_G(K);
            FloatVector sample_B(K);

            // points = IdentityMatrix(3); // TODO 
            
            for( int k = 0; k < K; k++ )
            {
                auto barycoords = get_random_barycentric_coordinates(2);
                assert( barycoords.isnonnegative() );
                assert( is_numerically_close( barycoords.sum(), 1. ) );
                
                auto P = M.getPointFromBarycentric( 2, t, barycoords );

                points.setrow( k, barycoords ); // TODO
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


        fstream fs( experimentfile( getbasename(__FILE__), "svg" ), std::fstream::out );
        fs << M.outputLinearSVG( interpol_red, interpol_green, interpol_blue, 0.000, "array", "none" );
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
        LOG << "Weights min/max: " << weights.front().second << space << weights.back().second << nl;
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
