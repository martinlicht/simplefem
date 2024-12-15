
#include <cstdio>

#include "../../basic.hpp"
#include "../../utility/math.hpp"

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for computing ball/sphere measures" << nl;

    for( int n = 0; n <= 10; n++ )
    {
        printf("Volume of %i-dimensional unit ball: %f\n", n, unitBallVolume(n) );
    }

    for( int n = 0; n <= 10; n++ )
    {
        printf("Area of %i-dimensional unit sphere: %f\n", n, unitSphereArea(n) );
    }

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
    return 0;
}
