#ifndef INCLUDEGUARD_UTILITY_MATH_HPP
#define INCLUDEGUARD_UTILITY_MATH_HPP


#include "../base/include.hpp"



// ============================================================================
// Bump functions
// ============================================================================

Float bumpfunction( Float x );

Float bumpfunction_dev( Float x );

Float bumpfunction_devdev( Float x );

Float bumpfunction_devdevdev( Float x );

Float blob( Float x );

Float blob_dev( Float x );

Float blob_devdev( Float x );

Float blob_devdevdev( Float x );

Float sinpy( Float x );

Float cospy( Float x );


// ============================================================================
// Cartesian and polar coordinates
// ============================================================================

void cartesian_to_polar_coordinates2D( const Float& x, const Float& y, Float& radius, Float& angle );

void polar_to_cartesian_coordinates2D( const Float& radius, const Float& angle, Float& x, Float& y );


// ============================================================================
// Volumes of unit ball, surface of unit sphere
// ============================================================================

Float unitBallVolume( int n );

Float unitSphereArea( int n );


#endif
