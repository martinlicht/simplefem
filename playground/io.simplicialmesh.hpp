#ifndef INCLUDEGUARD_IOSIMPLICIALMESH
#define INCLUDEGUARD_IOSIMPLICIALMESH


#include <string>
#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>


#include "../basic.hpp"
#include "coordinates.hpp"
#include "simplicialmesh.hpp"


/*******************
****  
****  Input / Output of simplicial meshes 
****  
*******************/


void writeSimplicialMesh( std::ostream&, const SimplicialMesh&, bool = false );

SimplicialMesh readSimplicialMesh( std::istream& );

void writeSimplicialMeshPath( const char*, const SimplicialMesh&, bool = false );

SimplicialMesh readSimplicialMeshPath( const char* );



#endif