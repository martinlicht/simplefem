
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <istream>
#include <ostream>
#include <fstream>



#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "mesh.simplicial1D.hpp"
#include "io.simplicial1D.hpp"
#include "coordinates.hpp"
#include "io.coordinates.hpp"




void writeMeshSimplicial1D( const char* filename, const MeshSimplicial1D& mesh, bool sugar )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    writeMeshSimplicial1D( myfile, mesh, sugar );
    myfile.close();
}

MeshSimplicial1D readMeshSimplicial1D( const char* filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    MeshSimplicial1D mesh = readMeshSimplicial1D( myfile );
    myfile.close();
    return mesh;
}



void writeMeshSimplicial1D( std::ostream& out, const MeshSimplicial1D& mesh, bool sugar )
{
    /* Preamble */
    if( sugar ) out << "Writing simplicial 1D Mesh..." << std::endl;
    
    if( sugar ) out << "number of edges: " << std::endl;;
    out << mesh.count_edges() << std::endl;
    
    if( sugar ) out << "number of vertices: " << std::endl;
    out << mesh.count_vertices() << std::endl;
    
    if( sugar ) out << "external dimension: " << std::endl;
    out << mesh.getcoordinates().getdimension() << std::endl;
    
    /* edge -> vertices */
    if( sugar ) out << "for each edge, the vertices: " << std::endl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_vertex( e, 0 )
            << space
            << mesh.get_edge_vertex( e, 1 ) << std::endl;
    }
    
    if( sugar ) out << "for each vertex, the first parent edge: " << std::endl;
    for( int v = 0; v < mesh.count_vertices(); v++ ) {
        if( sugar ) out << v << ": ";
        out << mesh.get_vertex_firstparent_edge( v )
            << std::endl;
    }
    
    /* edge -> next parent of vertex vertices */
    if( sugar ) out << "for each edge, the next neighbors: " << std::endl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_nextparent_of_vertex( e, 0 )
            << space
            << mesh.get_edge_nextparent_of_vertex( e, 1 )
            << std::endl;
    }
    
    assert( out.good() );
    
    writeCoordinates( out, mesh.getcoordinates(), sugar );
}




MeshSimplicial1D readMeshSimplicial1D( std::istream& in )
{
    int counter_edges, counter_vertices, dim;
    
    in >> counter_edges
       >> counter_vertices
       >> dim;
    
    int nullindex = MeshSimplicial1D::nullindex;
    
    /* edge -> vertices */
    
    std::vector< std::array<int,2> > edge_vertices( counter_edges, { nullindex, nullindex } );
    std::vector< int               > vertex_firstparent_edge( counter_vertices, nullindex );
    std::vector< std::array<int,2> > edge_nextparents_of_vertices( counter_edges, { nullindex, nullindex } );
    
    for( int e = 0; e < counter_edges; e++ )
        in >> edge_vertices[e][0] >> edge_vertices[e][1];
    
    for( int v = 0; v < counter_vertices; v++ )
        in >> vertex_firstparent_edge[v];
    
    for( int e = 0; e < counter_edges; e++ )
        in >> edge_nextparents_of_vertices[e][0] >> edge_nextparents_of_vertices[e][1];
    
    /* coordinates */
    
    Coordinates coords = readCoordinates( in );
    
    assert( in.good() );
    
    /* return */
    
    return MeshSimplicial1D( dim, coords, edge_vertices, edge_nextparents_of_vertices, vertex_firstparent_edge );    
}



        
        
        
        
        