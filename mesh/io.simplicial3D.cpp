
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
#include "mesh.simplicial3D.hpp"
#include "io.simplicial3D.hpp"
#include "coordinates.hpp"
#include "io.coordinates.hpp"




void writeMeshSimplicial3D( const char* filename, const MeshSimplicial3D& mesh, bool sugar )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    writeMeshSimplicial3D( myfile, mesh, sugar );
    myfile.close();
}

MeshSimplicial3D readMeshSimplicial3D( const char* filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    MeshSimplicial3D mesh = readMeshSimplicial3D( myfile );
    myfile.close();
    return mesh;
}



void writeMeshSimplicial3D( std::ostream& out, const MeshSimplicial3D& mesh, bool sugar )
{
    /* Preamble */
    if( sugar ) out << "Writing simplicial 3D Mesh..." << std::endl;
    
    if( sugar ) out << "number of tetrahedra: " << std::endl;;
    out << mesh.count_tetrahedra() << std::endl;
    
    if( sugar ) out << "number of faces: " << std::endl;;
    out << mesh.count_faces() << std::endl;
    
    if( sugar ) out << "number of edges: " << std::endl;;
    out << mesh.count_edges() << std::endl;
    
    if( sugar ) out << "number of vertices: " << std::endl;
    out << mesh.count_vertices() << std::endl;
    
    if( sugar ) out << "external dimension: " << std::endl;
    out << mesh.getcoordinates().getdimension() << std::endl;
    
    /* tetrahedron -> faces */
    
    if( sugar ) out << "for each tetrahedron, the faces: " << std::endl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_edge( t, 0 )
            << space
            << mesh.get_tetrahedron_edge( t, 1 )
            << space
            << mesh.get_tetrahedron_edge( t, 2 )
            << space
            << mesh.get_tetrahedron_edge( t, 3 ) << std::endl;
    }
    
    if( sugar ) out << "for each face, the first parent tetrahedron: " << std::endl;
    for( int f = 0; f < mesh.count_faces(); f++ ) {
        if( sugar ) out << f << ": ";
        out << mesh.get_face_firstparent_tetrahedron( f )
            << std::endl;
    }
    
    if( sugar ) out << "for each tetrahedron, the next neighbors: " << std::endl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_nextparent_of_face( t, 0 )
            << space
            << mesh.get_tetrahedron_nextparent_of_face( t, 1 )
            << space
            << mesh.get_tetrahedron_nextparent_of_face( t, 2 )
            << space
            << mesh.get_tetrahedron_nextparent_of_face( t, 3 )
            << std::endl;
    }
    
    /* tetrahedron -> edges */
    
    if( sugar ) out << "for each tetrahedron, the edges: " << std::endl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_edge( t, 0 )
            << space
            << mesh.get_tetrahedron_edge( t, 1 )
            << space
            << mesh.get_tetrahedron_edge( t, 2 )
            << space
            << mesh.get_tetrahedron_edge( t, 3 )
            << space
            << mesh.get_tetrahedron_edge( t, 4 )
            << space
            << mesh.get_tetrahedron_edge( t, 5 ) << std::endl;
    }
    
    if( sugar ) out << "for each edge, the first parent tetrahedron: " << std::endl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_firstparent_tetrahedron( e )
            << std::endl;
    }
    
    if( sugar ) out << "for each tetrahedron, the next neighbors: " << std::endl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_nextparent_of_edge( t, 0 )
            << space
            << mesh.get_tetrahedron_nextparent_of_edge( t, 1 )
            << space
            << mesh.get_tetrahedron_nextparent_of_edge( t, 2 )
            << space
            << mesh.get_tetrahedron_nextparent_of_edge( t, 3 )
            << space
            << mesh.get_tetrahedron_nextparent_of_edge( t, 4 )
            << space
            << mesh.get_tetrahedron_nextparent_of_edge( t, 5 )
            << std::endl;
    }
    
    /* tetrahedron -> vertices */
    
    if( sugar ) out << "for each tetrahedron, the vertices: " << std::endl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_vertex( t, 0 )
            << space
            << mesh.get_tetrahedron_vertex( t, 1 )
            << space
            << mesh.get_tetrahedron_vertex( t, 2 )
            << space
            << mesh.get_tetrahedron_vertex( t, 3 ) << std::endl;
    }
    
    if( sugar ) out << "for each vertex, the first parent tetrahedron: " << std::endl;
    for( int v = 0; v < mesh.count_vertices(); v++ ) {
        if( sugar ) out << v << ": ";
        out << mesh.get_vertex_firstparent_tetrahedron( v )
            << std::endl;
    }
    
    if( sugar ) out << "for each tetrahedron, the next neighbors: " << std::endl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_nextparent_of_vertex( t, 0 )
            << space
            << mesh.get_tetrahedron_nextparent_of_vertex( t, 1 )
            << space
            << mesh.get_tetrahedron_nextparent_of_vertex( t, 2 )
            << space
            << mesh.get_tetrahedron_nextparent_of_vertex( t, 3 )
            << std::endl;
    }
    
    assert( out.good() );
    
    /* face -> edges */
    
    if( sugar ) out << "for each face, the edges: " << std::endl;
    for( int f = 0; f < mesh.count_faces(); f++ ) {
        if( sugar ) out << f << ": ";
        out << mesh.get_face_edge( f, 0 )
            << space
            << mesh.get_face_edge( f, 1 )
            << space
            << mesh.get_face_edge( f, 2 ) << std::endl;
    }
    
    if( sugar ) out << "for each edge, the first parent face: " << std::endl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_firstparent_face( e )
            << std::endl;
    }
    
    if( sugar ) out << "for each face, the next neighbors: " << std::endl;
    for( int f = 0; f < mesh.count_faces(); f++ ) {
        if( sugar ) out << f << ": ";
        out << mesh.get_face_nextparent_of_edge( f, 0 )
            << space
            << mesh.get_face_nextparent_of_edge( f, 1 )
            << space
            << mesh.get_face_nextparent_of_edge( f, 2 )
            << std::endl;
    }
    
    /* face -> vertices */
    
    if( sugar ) out << "for each face, the vertices: " << std::endl;
    for( int f = 0; f < mesh.count_faces(); f++ ) {
        if( sugar ) out << f << ": ";
        out << mesh.get_face_vertex( f, 0 )
            << space
            << mesh.get_face_vertex( f, 1 )
            << space
            << mesh.get_face_vertex( f, 2 ) << std::endl;
    }
    
    if( sugar ) out << "for each vertex, the first parent face: " << std::endl;
    for( int v = 0; v < mesh.count_vertices(); v++ ) {
        if( sugar ) out << v << ": ";
        out << mesh.get_vertex_firstparent_face( v )
            << std::endl;
    }
    
    if( sugar ) out << "for each face, the next neighbors: " << std::endl;
    for( int f = 0; f < mesh.count_faces(); f++ ) {
        if( sugar ) out << f << ": ";
        out << mesh.get_face_nextparent_of_vertex( f, 0 )
            << space
            << mesh.get_face_nextparent_of_vertex( f, 1 )
            << space
            << mesh.get_face_nextparent_of_vertex( f, 2 )
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





MeshSimplicial3D readMeshSimplicial3D( std::istream& in )
{
    int counter_tetrahedra, counter_faces, counter_edges, counter_vertices, dim;
    
    in >> counter_tetrahedra
       >> counter_faces
       >> counter_edges
       >> counter_vertices
       >> dim;
    
    int nullindex = MeshSimplicial3D::nullindex;
    
    /* tetrahedron -> faces */
    
    std::vector< std::array<int,4> > tetrahedron_faces( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    std::vector< int               > face_firstparent_tetrahedron( counter_faces, nullindex );
    std::vector< std::array<int,4> > tetrahedron_nextparents_of_faces( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_faces[t][0]
           >> tetrahedron_faces[t][1]
           >> tetrahedron_faces[t][2]
           >> tetrahedron_faces[t][3];
    
    for( int f = 0; f < counter_faces; f++ )
        in >> face_firstparent_tetrahedron[f];
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_nextparents_of_faces[t][0] 
           >> tetrahedron_nextparents_of_faces[t][1] 
           >> tetrahedron_nextparents_of_faces[t][2] 
           >> tetrahedron_nextparents_of_faces[t][3];
    
    assert( in.good() );
    
    /* tetrahedron -> edges */
    
    std::vector< std::array<int,6> > tetrahedron_edges( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    std::vector< int               > edge_firstparent_tetrahedron( counter_edges, nullindex );
    std::vector< std::array<int,6> > tetrahedron_nextparents_of_edges( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_edges[t][0]
           >> tetrahedron_edges[t][1]
           >> tetrahedron_edges[t][2]
           >> tetrahedron_edges[t][3]
           >> tetrahedron_edges[t][4]
           >> tetrahedron_edges[t][5];
    
    for( int v = 0; v < counter_edges; v++ )
        in >> edge_firstparent_tetrahedron[v];
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_nextparents_of_edges[t][0] 
           >> tetrahedron_nextparents_of_edges[t][1] 
           >> tetrahedron_nextparents_of_edges[t][2] 
           >> tetrahedron_nextparents_of_edges[t][3] 
           >> tetrahedron_nextparents_of_edges[t][4] 
           >> tetrahedron_nextparents_of_edges[t][5];
    
    assert( in.good() );
    
    /* tetrahedron -> vertices */
    
    std::vector< std::array<int,4> > tetrahedron_vertices( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    std::vector< int               > vertex_firstparent_tetrahedron( counter_vertices, nullindex );
    std::vector< std::array<int,4> > tetrahedron_nextparents_of_vertices( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_vertices[t][0]
           >> tetrahedron_vertices[t][1]
           >> tetrahedron_vertices[t][2]
           >> tetrahedron_vertices[t][3];
    
    for( int v = 0; v < counter_vertices; v++ )
        in >> vertex_firstparent_tetrahedron[v];
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_nextparents_of_vertices[t][0] 
           >> tetrahedron_nextparents_of_vertices[t][1] 
           >> tetrahedron_nextparents_of_vertices[t][2] 
           >> tetrahedron_nextparents_of_vertices[t][3];
    
    assert( in.good() );
    
    /* face -> edges */
    
    std::vector< std::array<int,3> > face_edges( counter_faces, { nullindex, nullindex, nullindex } );
    std::vector< int               > edge_firstparent_face( counter_edges, nullindex );
    std::vector< std::array<int,3> > face_nextparents_of_edges( counter_faces, { nullindex, nullindex, nullindex } );
    
    for( int f = 0; f < counter_faces; f++ )
        in >> face_edges[f][0]
           >> face_edges[f][1]
           >> face_edges[f][2];
    
    for( int v = 0; v < counter_edges; v++ )
        in >> edge_firstparent_face[v];
    
    for( int f = 0; f < counter_faces; f++ )
        in >> face_nextparents_of_edges[f][0] 
           >> face_nextparents_of_edges[f][1] 
           >> face_nextparents_of_edges[f][2];
    
    assert( in.good() );
    
    /* face -> vertices */
    
    std::vector< std::array<int,3> > face_vertices( counter_faces, { nullindex, nullindex, nullindex } );
    std::vector< int               > vertex_firstparent_face( counter_vertices, nullindex );
    std::vector< std::array<int,3> > face_nextparents_of_vertices( counter_faces, { nullindex, nullindex, nullindex } );
    
    for( int f = 0; f < counter_faces; f++ )
        in >> face_vertices[f][0]
           >> face_vertices[f][1]
           >> face_vertices[f][2];
    
    for( int v = 0; v < counter_vertices; v++ )
        in >> vertex_firstparent_face[v];
    
    for( int f = 0; f < counter_faces; f++ )
        in >> face_nextparents_of_vertices[f][0] 
           >> face_nextparents_of_vertices[f][1] 
           >> face_nextparents_of_vertices[f][2];
    
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
    
    return MeshSimplicial3D( dim, coords, 
                             tetrahedron_faces, face_firstparent_tetrahedron, tetrahedron_nextparents_of_faces, 
                             tetrahedron_edges, edge_firstparent_tetrahedron, tetrahedron_nextparents_of_edges, 
                             tetrahedron_vertices, vertex_firstparent_tetrahedron, tetrahedron_nextparents_of_vertices, 
                             face_edges, edge_firstparent_face, face_nextparents_of_edges, 
                             face_vertices, vertex_firstparent_face, face_nextparents_of_vertices, 
                             edge_vertices, vertex_firstparent_edge, edge_nextparents_of_vertices
                           );    
}



