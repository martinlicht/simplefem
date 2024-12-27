
#include <cmath>
#include <cstdio>

#include "../../basic.hpp"
#include "../../utility/math.hpp"

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for computing ball/sphere measures" << nl;

    for( int n = 0; n <= 10; n++ )
    {
        LOGPRINTF("Area of %2i-dimensional unit sphere: %f\n", n, unitSphereArea(n) );
    }

    for( int n = 0; n <= 10; n++ )
    {
        LOGPRINTF("Volume of %2i-dimensional unit ball: %f\n", n, unitBallVolume(n) );
    }

    LOG << nl;

    #define M_PI Constants::pi

    // Reference values for volumes and areas up to dimension 10
    Float referenceAreas[10] = {
        2.0, 
        2.0 * M_PI, 
        4.0 * M_PI, 
        2.0 * M_PI * M_PI, 
        8.0 * M_PI * M_PI / 3.0, 
        //
        M_PI * M_PI * M_PI, 
        16.0 * M_PI * M_PI * M_PI / 15.0, 
        M_PI * M_PI * M_PI * M_PI / 3.0, 
        32.0 * pow(M_PI, 4) / 105.0, 
        pow(M_PI, 5) / 12.0
    };
    Float referenceVolumes[10] = {
        1.0, 
        2.0, 
        M_PI, 
        4. * M_PI / 3., 
        M_PI * M_PI / 2., 
        //
        8.0 * M_PI * M_PI / 15.0, 
        M_PI * M_PI * M_PI / 6.0, 
        16. * M_PI * M_PI * M_PI / 105., 
        M_PI * M_PI * M_PI * M_PI / 24.0, 
        32. * M_PI * M_PI * M_PI * M_PI / 945.0, 
        // M_PI * M_PI * M_PI * M_PI * M_PI / 120.0, 
    };



    // Loop through dimensions 0 to 10
    for (int n = 0; n < 10; n++) {

        Float computedArea = unitSphereArea(n);
        Float computedVolume = unitBallVolume(n);

        Float ra = absolute( computedArea   / referenceAreas[n]   );
        Float rv = absolute( computedVolume / referenceVolumes[n] );

        // Print results for verification
        std::cout << "Dimension " << n << ": Volume = " << computedVolume << " (Reference: " << referenceVolumes[n] << ")" << std::endl;
        std::cout << "Dimension " << n << ": Area = " << computedArea << " (Reference: " << referenceAreas[n] << ")" << std::endl;

        // Compare computed values with reference values
        Assert( is_numerically_one(ra,1e-10), n, computedArea,   referenceAreas[n]   );
        Assert( is_numerically_one(rv,1e-10), n, computedVolume, referenceVolumes[n] );

    }

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
    return 0;
}
