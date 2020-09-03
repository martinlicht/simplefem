#ifndef INCLUDEGUARD_MESH_IO_SIMPLICIALND_HPP
#define INCLUDEGUARD_MESH_IO_SIMPLICIALND_HPP


#include <istream>
#include <ostream>
#include <string>


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
