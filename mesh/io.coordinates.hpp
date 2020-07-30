#ifndef INCLUDEGUARD_MESH_IO_COORDINATES
#define INCLUDEGUARD_MESH_IO_COORDINATES


#include <istream>
#include <ostream>
#include <string>


#include "../basic.hpp"
#include "coordinates.hpp"


/*******************
****  
****  Input / Output of coordinate objects 
****  
*******************/


void writeCoordinates( std::ostream& out, const Coordinates& coordinates, bool sugar = false );

Coordinates readCoordinates( std::istream& in );

void writeCoordinates( const char* filename, const Coordinates& coordinates, bool sugar = false );

Coordinates readCoordinates( const char* filename );



#endif
