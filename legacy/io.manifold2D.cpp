
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>



#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "coordinates.hpp"
#include "io.coordinates.hpp"
#include "io.manifold2D.hpp"




void writeMeshManifold2D( const char* filename, const MeshManifold2D& mesh, bool sugar )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    writeMeshManifold2D( myfile, mesh, sugar );
    myfile.close();
}

MeshManifold2D readMeshManifold2D( const char* filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    MeshManifold2D mesh = readMeshManifold2D( myfile );
    myfile.close();
    return mesh;
}


void writeMeshManifold2D( std::ostream& out, const MeshManifold2D& mesh, bool sugar )
{
    /* Preamble */
    if( sugar ) out << "Mesh Manifold 2D" << std::endl;
    
    if( sugar ) out << "Triangles: \t"    << std::endl;
    out << mesh.count_triangles() << std::endl;
    if( sugar ) out << "Edges:     \t"    << std::endl;
    out << mesh.count_edges() << std::endl;
    if( sugar ) out << "Vertices:  \t"    << std::endl;
    out << mesh.count_vertices() << std::endl;
    
    
    /* data */
    for( int t = 0; t < mesh.count_triangles(); t++ ) {
        if( sugar ) out << t << ": ";
        for( int d = 0; d < 3; d++ )
            out << mesh.get_triangle_edges(t)[d] << "\t";
        out << std::endl;
    }
    
    for( int t = 0; t < mesh.count_triangles(); t++ ) {
        if( sugar ) out << t << ": ";
        for( int d = 0; d < 3; d++ )
            out << mesh.get_triangle_vertices(t)[d] << "\t";
        out << std::endl;
    }
    
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        for( int d = 0; d < 2; d++ )
            out << mesh.get_triangle_edge_parents(e)[d] << "\t";
        out << std::endl;
    }
    
    writeCoordinates( out, mesh.getcoordinates(), sugar );
    
}

MeshManifold2D readMeshManifold2D( std::istream& in )
{
    int counter_triangles;
    int counter_edges;
    int counter_vertices;
    in >> counter_triangles >> counter_edges >> counter_vertices;
    
    std::vector<std::array<int,3>> triangle_vertices( counter_triangles );
    std::vector<std::array<int,3>> triangle_edges( counter_triangles );
    std::vector<std::array<int,2>> edge_parents( counter_edges );
    
    for( int t = 0; t < counter_triangles; t++ )
        for( int d = 0; d < 3; d++ )
            in >> triangle_vertices[t][d];
    
    for( int t = 0; t < counter_triangles; t++ )
        for( int d = 0; d < 3; d++ )
            in >> triangle_edges[t][d];
    
    for( int e = 0; e < counter_edges; e++ )
        for( int d = 0; d < 2; d++ )
            in >> edge_parents[e][d];
    
    Coordinates coords = readCoordinates( in );
    
    return MeshManifold2D( 
              coords.getdimension(),
              coords,
              triangle_vertices,
              triangle_edges,
              edge_parents
           );
        
}


