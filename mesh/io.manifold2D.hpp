#ifndef INCLUDEGUARD_MESH_IO_MANIFOLD_2D
#define INCLUDEGUARD_MESH_IO_MANIFOLD_2D


#include <string>
#include <istream>
#include <ostream>


#include "../basic.hpp"
#include "mesh.manifold2D.hpp"


/*******************
****  
****  Input / Output of MeshManifold2D class
****  
*******************/


void writeMeshManifold2D( std::ostream& out, const MeshManifold2D& mm2d, bool sugar = false );

MeshManifold2D readMeshManifold2D( std::istream& in );

void writeMeshManifold2D( const char* filename, const MeshManifold2D& mm2d, bool sugar = false );

MeshManifold2D readMeshManifold2D( const char* filename );



#endif