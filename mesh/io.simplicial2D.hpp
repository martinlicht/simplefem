#ifndef INCLUDEGUARD_MESH_IO_SIMPLICIAL2D
#define INCLUDEGUARD_MESH_IO_SIMPLICIAL2D


#include <string>
#include <istream>
#include <ostream>


#include "../basic.hpp"
#include "mesh.simplicial2D.hpp"


/*******************
****  
****  Input / Output of coordinate objects 
****  
*******************/


void writeMeshSimplicial2D( std::ostream& out, const MeshSimplicial2D& mesh, bool sugar = false );

MeshSimplicial2D readMeshSimplicial2D( std::istream& in );

void writeMeshSimplicial2D( const char* filename, const MeshSimplicial2D& mesh, bool sugar = false );

MeshSimplicial2D readMeshSimplicial2D( const char* filename );



#endif