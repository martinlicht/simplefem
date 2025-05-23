
#include <cmath>

#include "../base/include.hpp"
#include "math.hpp"



/////////////////////////////////////////////////
//                                             //
//              BUMP FUNCTIONS                 //
//                                             //
/////////////////////////////////////////////////



Float bumpfunction( Float x )
{
    const Float delta = x * x - 1.;
    
    if( absolute( x ) < 0.99999999 ) {
    
        return std::exp( 1. / delta );
        
    } else {
    
        return 0;
        
    }
}

Float bumpfunction_dev( Float x )
{

    const Float delta = x * x - 1.;
    const Float delta_sq = delta * delta;
    
    if( absolute( x ) < 0.99999999 ) {
    
        return -2. * x * std::exp( 1. / delta ) / delta_sq;
        
    } else {
    
        return 0;
        
    }
}

Float bumpfunction_devdev( Float x )
{
    const Float delta = x * x - 1.;

    const Float x4 = x * x * x * x;
    
    const Float delta_sq = delta    * delta;
    const Float delta_p4 = delta_sq * delta_sq;
    
    if( absolute( x ) < 0.99999999 ) {
    
        return std::exp( 1. / delta ) * ( 6. * x4 - 2. ) / delta_p4;
        
    } else {
    
        return 0.;
        
    }
}



Float bumpfunction_devdevdev( Float x )
{
    const Float x2 = x * x;
    const Float x4 = x2 * x2;
    const Float x6 = x4 * x2;
    
    const Float delta = x2 - 1.;
    
    const Float delta_p2 = delta    * delta;
    const Float delta_p4 = delta_p2 * delta_p2;
    const Float delta_p6 = delta_p4 * delta_p2;
    
    if( absolute( x ) < 0.99999999 ) {
    
        return -4. * x * std::exp( 1. / delta ) * ( 6.*x6 + 3.*x4 - 10.*x2 + 3. ) / delta_p6;
        
    } else {
    
        return 0.;
        
    }
}



Float blob( Float x )
{
    const Float delta = x * x - 1.;
    
    if( absolute( x ) < 0.99999999 ) {
    
        return std::exp( 2. / delta );
        
    } else {
    
        return 0;
        
    }
}

Float blob_dev( Float x )
{

    const Float delta = x * x - 1.;
    const Float delta_sq = delta * delta;
    
    if( absolute( x ) < 0.99999999 ) {
    
        return -4. * x * std::exp( 2. / delta ) / delta_sq;
        
    } else {
    
        return 0;
        
    }
}

Float blob_devdev( Float x )
{

    const Float x2 = x * x;
    const Float x4 = x2 * x2;
    
    const Float delta = x2 - 1.;
    
    const Float delta_sq = delta    * delta;
    const Float delta_p4 = delta_sq * delta_sq;
    
    if( absolute( x ) < 0.99999999 ) {
    
        return 4. * std::exp( 2. / delta ) * ( 3.*x4 + 2.*x2 - 1. ) / delta_p4;
        
    } else {
    
        return 0.;
        
    }
}



Float blob_devdevdev( Float x )
{
    const Float x2 = x * x;
    const Float x3 = x2 * x;
    const Float x4 = x2 * x2;
    
    const Float delta = x2 - 1.;
    
    const Float delta_p2 = delta    * delta;
    const Float delta_p4 = delta_p2 * delta_p2;
    const Float delta_p6 = delta_p4 * delta_p2;
    
    if( absolute( x ) < 0.99999999 ) {
    
        return -16. * std::exp( 2. / delta ) * x3 * ( 3.*x4 + 6.*x2 - 5. ) / delta_p6;
        
    } else {
    
        return 0.;
        
    }
}



Float sinpy( Float x )
{
    return std::sin( Constants::pi * x );
}

Float cospy( Float x )
{
    return std::cos( Constants::pi * x );
}



/////////////////////////////////////////////////
//                                             //
//       CARTESIAN AND POLAR COORDINATES       //
//                                             //
/////////////////////////////////////////////////

void cartesian_to_polar_coordinates2D( const Float& x, const Float& y, Float& radius, Float& angle )
{
    radius = std::sqrt( x * x + y * y );
    angle  = std::atan2( y, x );
    
    if( angle < 0. ) angle = Constants::twopi + angle;
}

void polar_to_cartesian_coordinates2D( const Float& radius, const Float& angle, Float& x, Float& y )
{
    x = radius * std::cos( angle );
    y = radius * std::sin( angle );
}



// Function to compute the volume of the n-dimensional unit ball
Float unitBallVolume( int n )
{
    Assert( n >= 0 && "Dimension n must be non-negative" );
    
    if( n % 2 == 0 ) {
        // Even dimension: n = 2k
        int k = n / 2;
        return power_numerical( Constants::pi, k ) / factorial_numerical( k );
    } else {
        // Odd dimension: n = 2k+1
        int k = ( n - 1 ) / 2;
        return 2. * factorial_numerical( k ) * power_numerical( 4. * Constants::pi, k ) / factorial_numerical( 2 * k + 1 );
    }
}

// Function to compute the surface area of the n-dimensional unit sphere
Float unitSphereArea( int n )
{
    Assert( n >= 0 && "Dimension n must be non-negative" );
    return ( n + 1 ) * unitBallVolume( n + 1 );
}


