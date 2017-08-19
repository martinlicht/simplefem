
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
#include "mesh.simplicial2D.hpp"
#include "io.simplicial2D.hpp"
#include "coordinates.hpp"
#include "io.coordinates.hpp"




void writeMeshSimplicial2D( const char* filename, const MeshSimplicial2D& mesh, bool sugar )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    writeMeshSimplicial2D( myfile, mesh, sugar );
    myfile.close();
}

MeshSimplicial2D readMeshSimplicial2D( const char* filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    MeshSimplicial2D mesh = readMeshSimplicial2D( myfile );
    myfile.close();
    return mesh;
}



void writeMeshSimplicial2D( std::ostream& out, const MeshSimplicial2D& mesh, bool sugar )
{
    /* Preamble */
    if( sugar ) out << "Writing simplicial 2D Mesh..." << std::endl;
    
    if( sugar ) out << "number of triangles: " << std::endl;;
    out << mesh.count_triangles() << std::endl;
    
    if( sugar ) out << "number of edges: " << std::endl;;
    out << mesh.count_edges() << std::endl;
    
    if( sugar ) out << "number of vertices: " << std::endl;
    out << mesh.count_vertices() << std::endl;
    
    if( sugar ) out << "external dimension: " << std::endl;
    out << mesh.getcoordinates().getdimension() << std::endl;
    
    /* triangle -> edges */
    
    if( sugar ) out << "for each triangle, the edges: " << std::endl;
    for( int t = 0; t < mesh.count_triangles(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_triangle_edge( t, 0 )
            << space
            << mesh.get_triangle_edge( t, 1 )
            << space
            << mesh.get_triangle_edge( t, 2 ) << std::endl;
    }
    
    if( sugar ) out << "for each edge, the first parent triangle: " << std::endl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_firstparent_triangle( e )
            << std::endl;
    }
    
    if( sugar ) out << "for each triangle, the next neighbors: " << std::endl;
    for( int t = 0; t < mesh.count_triangles(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_triangle_nextparent_of_edge( t, 0 )
            << space
            << mesh.get_triangle_nextparent_of_edge( t, 1 )
            << space
            << mesh.get_triangle_nextparent_of_edge( t, 2 )
            << std::endl;
    }
    
    /* triangle -> vertices */
    
    if( sugar ) out << "for each triangle, the vertices: " << std::endl;
    for( int t = 0; t < mesh.count_triangles(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_triangle_vertex( t, 0 )
            << space
            << mesh.get_triangle_vertex( t, 1 )
            << space
            << mesh.get_triangle_vertex( t, 2 ) << std::endl;
    }
    
    if( sugar ) out << "for each vertex, the first parent triangle: " << std::endl;
    for( int v = 0; v < mesh.count_vertices(); v++ ) {
        if( sugar ) out << v << ": ";
        out << mesh.get_vertex_firstparent_triangle( v )
            << std::endl;
    }
    
    if( sugar ) out << "for each triangle, the next neighbors: " << std::endl;
    for( int t = 0; t < mesh.count_triangles(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_triangle_nextparent_of_vertex( t, 0 )
            << space
            << mesh.get_triangle_nextparent_of_vertex( t, 1 )
            << space
            << mesh.get_triangle_nextparent_of_vertex( t, 2 )
            << std::endl;
    }
    
    assert( out.good() );
    
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





MeshSimplicial2D readMeshSimplicial2D( std::istream& in )
{
    int counter_triangles, counter_edges, counter_vertices, dim;
    
    in >> counter_triangles
       >> counter_edges
       >> counter_vertices
       >> dim;
    
    int nullindex = MeshSimplicial2D::nullindex;
    
    /* triangle -> edges */
    
    std::vector< std::array<int,3> > triangle_edges( counter_triangles, { nullindex, nullindex, nullindex } );
    std::vector< int               > edge_firstparent_triangle( counter_edges, nullindex );
    std::vector< std::array<int,3> > triangle_nextparents_of_edges( counter_triangles, { nullindex, nullindex, nullindex } );
    
    for( int t = 0; t < counter_triangles; t++ )
        in >> triangle_edges[t][0]
           >> triangle_edges[t][1]
           >> triangle_edges[t][2];
    
    for( int v = 0; v < counter_edges; v++ )
        in >> edge_firstparent_triangle[v];
    
    for( int t = 0; t < counter_triangles; t++ )
        in >> triangle_nextparents_of_edges[t][0] 
           >> triangle_nextparents_of_edges[t][1] 
           >> triangle_nextparents_of_edges[t][2];
    
    assert( in.good() );
    
    /* triangle -> vertices */
    
    std::vector< std::array<int,3> > triangle_vertices( counter_triangles, { nullindex, nullindex, nullindex } );
    std::vector< int               > vertex_firstparent_triangle( counter_vertices, nullindex );
    std::vector< std::array<int,3> > triangle_nextparents_of_vertices( counter_triangles, { nullindex, nullindex, nullindex } );
    
    for( int t = 0; t < counter_triangles; t++ )
        in >> triangle_vertices[t][0]
           >> triangle_vertices[t][1]
           >> triangle_vertices[t][2];
    
    for( int v = 0; v < counter_vertices; v++ )
        in >> vertex_firstparent_triangle[v];
    
    for( int t = 0; t < counter_triangles; t++ )
        in >> triangle_nextparents_of_vertices[t][0] 
           >> triangle_nextparents_of_vertices[t][1] 
           >> triangle_nextparents_of_vertices[t][2];
    
    assert( in.good() );
    
    /* edges -> vertices */
    
    std::vector< std::array<int,2> > edge_vertices( counter_edges, { nullindex, nullindex } );
    std::vector< int               > vertex_firstparent_edge( counter_vertices, nullindex );
    std::vector< std::array<int,2> > edge_nextparents_of_vertices( counter_edges, { nullindex, nullindex } );
    
    for( int e = 0; e < counter_edges; e++ )
        in >> edge_vertices[e][0] >> edge_vertices[e][1];
    
    for( int v = 0; v < counter_vertices; v++ )
        in >> vertex_firstparent_edge[v];
    
    for( int e = 0; e < counter_edges; e++ )
        in >> edge_nextparents_of_vertices[e][0] >> edge_nextparents_of_vertices[e][1];
    
    assert( in.good() );
    
    /* coordinates */
    
    Coordinates coords = readCoordinates( in );
    
    assert( in.good() );
    
    /* return */
    
    return MeshSimplicial2D( dim, coords, 
                             triangle_edges, edge_firstparent_triangle, triangle_nextparents_of_edges, 
                             triangle_vertices, vertex_firstparent_triangle, triangle_nextparents_of_vertices, 
                             edge_vertices, vertex_firstparent_edge, edge_nextparents_of_vertices
                           );    
}



