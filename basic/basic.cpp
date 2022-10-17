
#include <chrono>
#include <vector>

#include "basic.hpp"
#include "constants.hpp"

template class std::vector<int>;
template class std::vector<Float>;


// Since all literals throughout are double unless marked otherwise 
// we enforce that `Float` is at least enough to store double.
// Any of those should do:
// 
// static_assert( Float(std::numeric_limits<double>::max()) == std::numeric_limits<double>::max(), "Float must be at least double" );
static_assert( sizeof(Float) >= sizeof(double), "Float must be at least double" );


#include <cstdarg>

std::string printf_into_string( const char* formatstring, ... )
{
    
    va_list args;
    
    va_start( args, formatstring );
    std::size_t length = std::vsnprintf(nullptr, 0, formatstring, args ) + 1;
    va_end( args );
    
    char* c_str = new char[length];
    
    va_start( args, formatstring );
    std::vsnprintf( c_str, length, formatstring, args );
    va_end( args );
    
    std::string ret( c_str );
    delete[] c_str;
    return ret;
}





// TODO: Move time to cpp
static_assert( std::is_integral< decltype( std::chrono::time_point_cast< std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count() ) >::value , "Time measurement must be integral" );

timestamp gettimestamp()
{
    
    static timestamp start_timestamp = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    timestamp now = std::chrono::time_point_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now() ).time_since_epoch().count();
    
    Assert( now >= start_timestamp );
    
    return now - start_timestamp;
}



// TODO: move to utility 

std::string timestamp2measurement( const timestamp& t )
{
    return std::to_string( static_cast<uintmax_t>(t) ) + "ms";
}

// std::string measurementnow( const timestamp& t ) // TODO Remove this line 
std::string measurementnow()
{
    return timestamp2measurement( gettimestamp() );
}



std::string timestamp2digitalcode( const timestamp& t )
{
    const int numdigits = 10;
    const int fulllength = numdigits+1;
    char digits[fulllength];
    snprintf( digits, fulllength, "%*ju", numdigits, (uintmax_t)t );
    for( int i = 0; i < numdigits; i++ ) if( digits[i] == ' ' ) digits[i] = '_';
    return std::string(digits);
}

std::string digitalcodenow()
{
    return timestamp2digitalcode( gettimestamp() );
}



std::string protocolprefixnow()
{
    // static const std::string foo = std::string("\e[36m[");
    // static const std::string bar = std::string("]\e[39m\t");
    static const std::string foo = std::string("[");
    static const std::string bar = std::string("]\t");
    return foo + digitalcodenow() + bar;
}










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

