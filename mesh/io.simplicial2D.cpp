
#include <algorithm>
#include <fstream>
#include <istream>
#include <ostream>
#include <map>
#include <string>
#include <vector>
#include <utility>



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
    if( sugar ) out << "Writing simplicial 2D Mesh..." << nl;
    
    if( sugar ) out << "number of triangles: " << nl;;
    out << mesh.count_triangles() << nl;
    
    if( sugar ) out << "number of edges: " << nl;;
    out << mesh.count_edges() << nl;
    
    if( sugar ) out << "number of vertices: " << nl;
    out << mesh.count_vertices() << nl;
    
    if( sugar ) out << "external dimension: " << nl;
    out << mesh.getcoordinates().getdimension() << nl;
    
    /* triangle -> edges */
    
    if( sugar ) out << "for each triangle, the edges: " << nl;
    for( int t = 0; t < mesh.count_triangles(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_triangle_edge( t, 0 )
            << space
            << mesh.get_triangle_edge( t, 1 )
            << space
            << mesh.get_triangle_edge( t, 2 ) << nl;
    }
    
    if( sugar ) out << "for each edge, the first parent triangle: " << nl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_firstparent_triangle( e )
            << nl;
    }
    
    if( sugar ) out << "for each triangle, the next neighbors: " << nl;
    for( int t = 0; t < mesh.count_triangles(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_triangle_nextparent_of_edge( t, 0 )
            << space
            << mesh.get_triangle_nextparent_of_edge( t, 1 )
            << space
            << mesh.get_triangle_nextparent_of_edge( t, 2 )
            << nl;
    }
    
    /* triangle -> vertices */
    
    if( sugar ) out << "for each triangle, the vertices: " << nl;
    for( int t = 0; t < mesh.count_triangles(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_triangle_vertex( t, 0 )
            << space
            << mesh.get_triangle_vertex( t, 1 )
            << space
            << mesh.get_triangle_vertex( t, 2 ) << nl;
    }
    
    if( sugar ) out << "for each vertex, the first parent triangle: " << nl;
    for( int v = 0; v < mesh.count_vertices(); v++ ) {
        if( sugar ) out << v << ": ";
        out << mesh.get_vertex_firstparent_triangle( v )
            << nl;
    }
    
    if( sugar ) out << "for each triangle, the next neighbors: " << nl;
    for( int t = 0; t < mesh.count_triangles(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_triangle_nextparent_of_vertex( t, 0 )
            << space
            << mesh.get_triangle_nextparent_of_vertex( t, 1 )
            << space
            << mesh.get_triangle_nextparent_of_vertex( t, 2 )
            << nl;
    }
    
    assert( out.good() );
    
    /* edge -> vertices */
    
    if( sugar ) out << "for each edge, the vertices: " << nl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_vertex( e, 0 )
            << space
            << mesh.get_edge_vertex( e, 1 ) << nl;
    }
    
    if( sugar ) out << "for each vertex, the first parent edge: " << nl;
    for( int v = 0; v < mesh.count_vertices(); v++ ) {
        if( sugar ) out << v << ": ";
        out << mesh.get_vertex_firstparent_edge( v )
            << nl;
    }
    
    if( sugar ) out << "for each edge, the next neighbors: " << nl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_nextparent_of_vertex( e, 0 )
            << space
            << mesh.get_edge_nextparent_of_vertex( e, 1 )
            << nl;
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



