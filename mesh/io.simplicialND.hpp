#ifndef INCLUDEGUARD_MESH_IO_SIMPLICIALND
#define INCLUDEGUARD_MESH_IO_SIMPLICIALND


#include <string>
#include <istream>
#include <ostream>


#include "../basic.hpp"
#include "mesh.simplicialND.hpp"


/*******************
****  
****  Input / Output of coordinate objects 
****  
*******************/


void writeMeshSimplicialND( std::ostream& out, const MeshSimplicialND& mesh, bool sugar = false );

MeshSimplicialND readMeshSimplicialND( std::istream& in );

void writeMeshSimplicialND( const char* filename, const MeshSimplicialND& mesh, bool sugar = false );

MeshSimplicialND readMeshSimplicialND( const char* filename );



#endif