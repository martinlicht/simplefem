#ifndef INCLUDEGUARD_IOCOORDINATES
#define INCLUDEGUARD_IOCOORDINATES


#include <string>
#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>


#include "../basic.hpp"
#include "coordinates.hpp"


/*******************
****  
****  Input / Output of simplicial meshes 
****  
*******************/


void writeCoordinates( std::ostream&, const Coordinates&, bool = false );

Coordinates readCoordinates( std::istream& );

void writeCoordinatesPath( const char*, const Coordinates&, bool = false );

Coordinates readCoordinatesPath( const char* );



#endif