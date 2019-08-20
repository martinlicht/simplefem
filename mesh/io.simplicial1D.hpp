#ifndef INCLUDEGUARD_MESH_IO_SIMPLICIAL1D
#define INCLUDEGUARD_MESH_IO_SIMPLICIAL1D


#include <string>
#include <istream>
#include <ostream>


#include "../basic.hpp"
#include "mesh.simplicial1D.hpp"


/*******************
****  
****  Input / Output of coordinate objects 
****  
*******************/


void writeMeshSimplicial1D( std::ostream& out, const MeshSimplicial1D& mesh, bool sugar = false );

MeshSimplicial1D readMeshSimplicial1D( std::istream& in );

void writeMeshSimplicial1D( const char* filename, const MeshSimplicial1D& mesh, bool sugar = false );

MeshSimplicial1D readMeshSimplicial1D( const char* filename );



#endif