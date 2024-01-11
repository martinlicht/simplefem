
#include "../basic/basic.hpp"
#include "math.hpp"



/////////////////////////////////////////////////
//                                             //
//              BUMP FUNCTIONS                 //
//                                             //
/////////////////////////////////////////////////





Float bumpfunction( Float x )
{
    Float delta = x*x - 1.;

    if( absolute(x) < 0.99999999 ) {

        return std::exp( 1. / delta );

    } else {

        return 0;
        
    }
}

Float bumpfunction_dev( Float x )
{
    
    Float delta = x*x - 1.;
    Float delta_sq = delta*delta;

    if( absolute(x) < 0.99999999 ) {
        
        return -2. * x * std::exp( 1. / delta ) / delta_sq;
        
    } else {
        
        return 0;
        
    }
}

Float bumpfunction_devdev( Float x )
{
    
    
//     Float t1 = std::exp( -1. / ( 1. - x*x ) );
//     Float t2 = std::exp( 1 - x*x );
// 
//     if( x*x < 1 )
//         return
//             t1 * ( 4*x*x*pow(t2,-4.) - 2*pow(t2,-2.) - 8*x*x*pow(t2,-3.) );
//     else
//         return
//             0.;
                            

    
    
    Float delta = x*x - 1.;
    
    Float delta_sq = delta    * delta;
    Float delta_p4 = delta_sq * delta_sq;

    if( absolute(x) < 0.99999999 ) {
        
        return std::exp( 1. / delta ) * (  6.*x*x*x*x - 2. ) / delta_p4;
        
    } else {
        
        return 0.;
        
    }
}




Float bumpfunction_devdevdev( Float x )
{
    // TODO: compute the correct value     

    Float x2 = x*x;
    Float x4 = x2 * x2;
    Float x6 = x4 * x2;

    Float delta = x2 - 1.;
    
    Float delta_p2 = delta    * delta;
    Float delta_p4 = delta_p2 * delta_p2;
    Float delta_p6 = delta_p4 * delta_p2;

    if( absolute(x) < 0.99999999 ) {
        
        return -4. * x * std::exp( 1. / delta ) * (  6.*x6 + 3.*x4 - 10.*x2 + 3. ) / delta_p6;
        
    } else {
        
        return 0.;
        
    }
}





Float blob( Float x )
{
    Float delta = x*x - 1.;

    if( absolute(x) < 0.99999999 ) {

        return std::exp( 2. / delta );

    } else {

        return 0;
        
    }
}

Float blob_dev( Float x )
{
    
    Float delta = x*x - 1.;
    Float delta_sq = delta*delta;

    if( absolute(x) < 0.99999999 ) {
        
        return -4. * x * std::exp( 2. / delta ) / delta_sq;
        
    } else {
        
        return 0;
        
    }
}

Float blob_devdev( Float x )
{

    Float x2 = x*x;
    Float x4 = x2 * x2;

    Float delta = x2 - 1.;
    
    Float delta_sq = delta    * delta;
    Float delta_p4 = delta_sq * delta_sq;

    if( absolute(x) < 0.99999999 ) {
        
        return 4. * std::exp( 2. / delta ) * (  3.*x4 + 2.*x2 - 1. ) / delta_p4;
        
    } else {
        
        return 0.;
        
    }
}




Float blob_devdevdev( Float x )
{
    // TODO: compute the correct value     

    Float x2 = x*x;
    Float x3 = x2*x;
    Float x4 = x2 * x2;
    
    Float delta = x2 - 1.;
    
    Float delta_p2 = delta    * delta;
    Float delta_p4 = delta_p2 * delta_p2;
    Float delta_p6 = delta_p4 * delta_p2;

    if( absolute(x) < 0.99999999 ) {
        
        return -16. * std::exp( 2. / delta ) * x3 * (  3.*x4 + 6.*x2 - 5. ) / delta_p6;
        
    } else {
        
        return 0.;
        
    }
}





Float sinpy( Float x )
{
    return sin( Constants::pi * x );
}

Float cospy( Float x )
{
    return cos( Constants::pi * x );
}




/////////////////////////////////////////////////
//                                             //
//       CARTESIAN AND POLAR COORDINATES       //
//                                             //
/////////////////////////////////////////////////

void cartesian_to_polar_coordinates2D( const Float& x, const Float& y, Float& radius, Float& angle )
{
    radius = std::sqrt( x*x + y*y );
    angle  = std::atan2( y, x );
    if( angle < 0. ) angle = Constants::twopi + angle;
}

void polar_to_cartesian_coordinates2D( const Float& radius, const Float& angle, Float& x, Float& y )
{
    x = radius * std::cos( angle );
    y = radius * std::sin( angle );
}

