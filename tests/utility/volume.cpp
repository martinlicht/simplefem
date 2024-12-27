
#include <cmath>
#include <cstdio>

#include "../../basic.hpp"
#include "../../utility/math.hpp"

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for computing ball/sphere measures" << nl;

    const Float pi = Constants::pi;

    // Reference values for areas and volumes up to dimension 9
    Float referenceAreas[10] = {
        2.0, 
        2.0 * pi, 
        4.0 * pi, 
        2.0 * pi * pi, 
        8.0 * pi * pi / 3.0, 
        //
        pi * pi * pi, 
        16.0 * pi * pi * pi / 15.0, 
        pi * pi * pi * pi / 3.0, 
        32.0 * pow(pi, 4) / 105.0, 
        pow(pi, 5) / 12.0
    };

    Float referenceVolumes[10] = {
        1.0, 
        2.0, 
        pi, 
        4. * pi / 3., 
        pi * pi / 2., 
        //
        8.0 * pi * pi / 15.0, 
        pi * pi * pi / 6.0, 
        16. * pi * pi * pi / 105., 
        pi * pi * pi * pi / 24.0, 
        32. * pi * pi * pi * pi / 945.0, 
        // pi * pi * pi * pi * pi / 120.0, 
    };



    for( int n = 0; n <= 9; n++ )
    {
        Float computedArea = unitSphereArea(n);
        LOGPRINTF("Dimension %2i-dimensional area sphere: %.10f / %.10f\n", n, unitSphereArea(n), referenceAreas[n] );
        Float ra = absolute( computedArea / referenceAreas[n] );
        Assert( is_numerically_one(ra,1e-10), n, computedArea,   referenceAreas[n]   );
    }

    for( int n = 0; n <= 9; n++ )
    {
        Float computedVolume = unitBallVolume(n);
        LOGPRINTF("Dimension %2i-dimensional volume ball: %.10f / %.10f\n", n, unitBallVolume(n), referenceVolumes[n] );
        Float rv = absolute( computedVolume / referenceVolumes[n] );
        Assert( is_numerically_one(rv,1e-10), n, computedVolume, referenceVolumes[n] );
    }

    LOG << nl;

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
    return 0;
}
