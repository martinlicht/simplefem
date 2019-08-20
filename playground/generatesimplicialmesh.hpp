#ifndef INCLUDEGUARD_GENERATESIMPLICIALMESH
#define INCLUDEGUARD_GENERATESIMPLICIALMESH


#include <string>
#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>


#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../operators/floatvector.hpp"
#include "coordinates.hpp"
#include "simplicialmesh.hpp"



/*******************
****  
****  Generate example meshes 
****  
****  - Unit cube, any dimension, simple triangulation
****  - Unit cube, 2D + 3D, resolution can be controlled 
****  - Load mesh from stream 
****  
*******************/




SimplicialMesh UnitCubeTriangulation( int innerdim, int outerdim );

SimplicialMesh UnitCubeTriangulation2D( int xstep, int ystep );

SimplicialMesh UnitCubeTriangulation3D( int xstep, int ystep, int zstep );

void generateMeshFromStream( std::istream& in );



#endif