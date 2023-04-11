

/**/
#include <algorithm>
#include <vector>

#include "../../basic.hpp"
#include "../../utility/files.hpp"
#include "../../utility/pixelimage.hpp"
#include "../../utility/random.hpp"
#include "../../dense/qr.factorization.hpp"
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

    // std::string name = "aurora.jpeg";
    std::string name = "lena_color.tiff";
    // std::string name = "testbild.jpg";

    // PixelImage pim = readPixelImage("lena_color.tiff");
    PixelImage pim = readPixelImage( name );

    auto red   = pim.get_interpolated_red();
    auto green = pim.get_interpolated_green();
    auto blue  = pim.get_interpolated_blue();
    
    M.getcoordinates().scale( { (Float)pim.getwidth(), (Float)pim.getheight() } );

    int l_min =  0;
    int l_max =  15;
    
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

            for( int k = 0; k < K; k++ )
            {
                auto barycoords = get_random_barycentric_coordinates(2);
                auto P = M.getPointFromBarycentric( 2, t, barycoords );

                points.setrow( k, barycoords );
                sample_R[k] = red(   P[0], P[1] );
                sample_G[k] = green( P[0], P[1] );
                sample_B[k] = blue(  P[0], P[1] );
            }

            auto local_red   = SolveOverconstrained( points, sample_R );
            auto local_green = SolveOverconstrained( points, sample_G );
            auto local_blue  = SolveOverconstrained( points, sample_B );

            for( int c = 0; c <= 2; c++ ) {
                interpol_red[  3*t+c] = sample_R[c];
                interpol_green[3*t+c] = sample_G[c];
                interpol_blue[ 3*t+c] = sample_B[c];
            }
        }


        // FloatVector redvec( num_tets, 128 ), greenvec( num_tets, 240 ), bluevec( num_tets, 38 );
        // redvec.random_within_range(0.,255.); greenvec.random_within_range(0.,255.); blue.random_within_range(0.,255.);

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
                samples_R[k] = red(   P[0], P[1] ) - barycoords[0] * interpol_red[3*t+0]   - barycoords[1] * interpol_red[3*t+1]   - barycoords[2] * interpol_red[3*t+2];
                samples_G[k] = green( P[0], P[1] ) - barycoords[0] * interpol_green[3*t+0] - barycoords[1] * interpol_green[3*t+1] - barycoords[2] * interpol_green[3*t+2];
                samples_B[k] = blue(  P[0], P[1] ) - barycoords[0] * interpol_blue[3*t+0]  - barycoords[1] * interpol_blue[3*t+1]  - barycoords[2] * interpol_blue[3*t+2];
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

    LOG << "Finished Unit Test" << nl;

    return 0;
}
