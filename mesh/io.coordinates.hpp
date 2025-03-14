#ifndef INCLUDEGUARD_MESH_IO_COORDINATES_HPP
#define INCLUDEGUARD_MESH_IO_COORDINATES_HPP


#include <istream>
#include <ostream>


#include "../base/include.hpp"
#include "coordinates.hpp"


/*******************
****  
****  Input / Output of coordinate objects 
****  
*******************/


void writeCoordinates( std::ostream& out, const Coordinates& coords, bool sugar = false );

Coordinates readCoordinates( std::istream& in );

void writeCoordinates( const char* filename, const Coordinates& coords, bool sugar = false );

Coordinates readCoordinates( const char* filename );



#endif
