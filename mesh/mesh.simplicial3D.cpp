
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
#include "mesh.hpp"

#include "mesh.simplicial3D.hpp"




MeshSimplicial3D::MeshSimplicial3D( int outerdim )
:
    Mesh( 3, outerdim ),
    
    counter_tetrahedra(0),
    counter_faces(0),
    counter_edges(0),
    counter_vertices(0),
    
    data_tetrahedron_faces(0),
    data_face_firstparent_tetrahedron(0),
    data_tetrahedron_nextparents_of_faces(0),
    
    data_tetrahedron_edges(0),
    data_edge_firstparent_tetrahedron(0),
    data_tetrahedron_nextparents_of_edges(0),
    
    data_tetrahedron_vertices(0),
    data_vertex_firstparent_tetrahedron(0),
    data_tetrahedron_nextparents_of_vertices(0),
    
    data_face_edges(0),
    data_edge_firstparent_face(0),
    data_face_nextparents_of_edges(0),
    
    data_face_vertices(0),
    data_vertex_firstparent_face(0),
    data_face_nextparents_of_vertices(0),
    
    data_edge_vertices(0),
    data_vertex_firstparent_edge(0),
    data_edge_nextparents_of_vertices(0)
{
    check();
}


MeshSimplicial3D::MeshSimplicial3D( 
    int outerdim,
    const Coordinates& coords,
    const std::vector<std::array<int,4>> tetrahedron_vertices
)
:
    Mesh( 3, outerdim ),
    
    counter_tetrahedra( tetrahedron_vertices.size() ),
    counter_faces( 0 ),
    counter_edges( 0 ),
    counter_vertices( 0 ),
    
    data_tetrahedron_faces( tetrahedron_vertices.size(), { nullindex, nullindex, nullindex, nullindex } ),
    data_face_firstparent_tetrahedron( 0 ),
    data_tetrahedron_nextparents_of_faces( tetrahedron_vertices.size(), { nullindex, nullindex, nullindex, nullindex } ),
    
    data_tetrahedron_edges( tetrahedron_vertices.size(), { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } ),
    data_edge_firstparent_tetrahedron( 0 ),
    data_tetrahedron_nextparents_of_edges( tetrahedron_vertices.size(), { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } ),
    
    data_tetrahedron_vertices( tetrahedron_vertices ),
    data_vertex_firstparent_tetrahedron( 0 ),
    data_tetrahedron_nextparents_of_vertices( tetrahedron_vertices.size(), { nullindex, nullindex, nullindex, nullindex } ),
    
    data_face_edges( 0, { nullindex, nullindex, nullindex } ),
    data_edge_firstparent_face( 0 ),
    data_face_nextparents_of_edges( tetrahedron_vertices.size(), { nullindex, nullindex, nullindex } ),
    
    data_face_vertices( 0, { nullindex, nullindex, nullindex } ),
    data_vertex_firstparent_face( 0 ),
    data_face_nextparents_of_vertices( tetrahedron_vertices.size(), { nullindex, nullindex, nullindex } ),
    
    data_edge_vertices( 0 ),
    data_vertex_firstparent_edge( 0 ),
    data_edge_nextparents_of_vertices( 0 )
{
    
    getcoordinates() = coords;
    
    /* 1. create all faces */
    
    data_face_vertices.resize( counter_tetrahedra * 4 );
    for( int t  = 0; t  < counter_tetrahedra;  t++ )
    {
      data_face_vertices[ 0 * counter_tetrahedra + t ] = { data_tetrahedron_vertices[t][0], data_tetrahedron_vertices[t][1], data_tetrahedron_vertices[t][2] };
      data_face_vertices[ 1 * counter_tetrahedra + t ] = { data_tetrahedron_vertices[t][0], data_tetrahedron_vertices[t][1], data_tetrahedron_vertices[t][3] };
      data_face_vertices[ 2 * counter_tetrahedra + t ] = { data_tetrahedron_vertices[t][0], data_tetrahedron_vertices[t][2], data_tetrahedron_vertices[t][3] };
      data_face_vertices[ 3 * counter_tetrahedra + t ] = { data_tetrahedron_vertices[t][1], data_tetrahedron_vertices[t][2], data_tetrahedron_vertices[t][3] };
    }
    
    {

    /* Due to an obscure bug in the C++ standard itself (!), we cannot use the standard sort function here */
    /* For that reason, a hand-written sort is used. No complexity estimates are guaranteed. */
    
//       std::sort( data_face_vertices.begin(), data_face_vertices.end() ); 
      sorthack( data_face_vertices );
      auto it = std::unique( data_face_vertices.begin(), data_face_vertices.end() );
      data_face_vertices.resize( it - data_face_vertices.begin() );
    }
    
    counter_faces = data_face_vertices.size();
    
    /* 2. create all edges */
    
    data_edge_vertices.resize( counter_tetrahedra * 6 );
    for( int t  = 0; t  < counter_tetrahedra;  t++ )
    {
      data_edge_vertices[ 0 * counter_tetrahedra + t ] = { data_tetrahedron_vertices[t][0], data_tetrahedron_vertices[t][1] };
      data_edge_vertices[ 1 * counter_tetrahedra + t ] = { data_tetrahedron_vertices[t][0], data_tetrahedron_vertices[t][2] };
      data_edge_vertices[ 2 * counter_tetrahedra + t ] = { data_tetrahedron_vertices[t][0], data_tetrahedron_vertices[t][3] };
      data_edge_vertices[ 3 * counter_tetrahedra + t ] = { data_tetrahedron_vertices[t][1], data_tetrahedron_vertices[t][2] };
      data_edge_vertices[ 4 * counter_tetrahedra + t ] = { data_tetrahedron_vertices[t][1], data_tetrahedron_vertices[t][3] };
      data_edge_vertices[ 5 * counter_tetrahedra + t ] = { data_tetrahedron_vertices[t][2], data_tetrahedron_vertices[t][3] };
    }
    
    {

    /* Due to an obscure bug in the C++ standard itself (!), we cannot use the standard sort function here */
    /* For that reason, a hand-written sort is used. No complexity estimates are guaranteed. */
    
//       std::sort( data_edge_vertices.begin(), data_edge_vertices.end() ); 
      sorthack( data_edge_vertices );
      auto it = std::unique( data_edge_vertices.begin(), data_edge_vertices.end() );
      data_edge_vertices.resize( it - data_edge_vertices.begin() );
    }
    
    counter_edges = data_edge_vertices.size();
    
    /* 3. Count vertices */
    
    counter_vertices = 0;
    for( const auto& quartett : data_tetrahedron_vertices )
    for( const int& vertex : quartett )
      counter_vertices = counter_vertices < vertex ? vertex : counter_vertices; 
    counter_vertices += 1;
    
    
    /* 4. Allocate the remaining memory */
    
    data_tetrahedron_faces.resize( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_edges.resize( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    data_face_edges.resize( counter_faces, { nullindex, nullindex, nullindex } );
    
    
    data_tetrahedron_nextparents_of_faces.resize( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_nextparents_of_edges.resize( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_nextparents_of_vertices.resize( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    data_face_nextparents_of_edges.resize( counter_faces, { nullindex, nullindex, nullindex } );
    data_face_nextparents_of_vertices.resize( counter_faces, { nullindex, nullindex, nullindex } );
    
    data_edge_nextparents_of_vertices.resize( counter_edges, { nullindex, nullindex } );
    
    
    data_face_firstparent_tetrahedron.resize( counter_faces, nullindex );
    
    data_edge_firstparent_tetrahedron.resize( counter_edges, nullindex );
    data_edge_firstparent_face.resize( counter_edges, nullindex );
    
    data_vertex_firstparent_tetrahedron.resize( counter_vertices, nullindex );
    data_vertex_firstparent_face.resize( counter_vertices, nullindex );
    data_vertex_firstparent_edge.resize( counter_vertices, nullindex );
    
    
    /* 5. For each vertex, set the first parent tetrahedron and the neighboring parent tetrahedra */
    
    for( int t =  0; t  < counter_tetrahedra; t++  )
    for( int vi = 0; vi <                  4; vi++ )
    {
      int v = data_tetrahedron_vertices[t][vi];
      
      assert( 0 <= v && v < counter_vertices );
      
      if( data_vertex_firstparent_tetrahedron[v] == nullindex ) {
        
        data_vertex_firstparent_tetrahedron[v] = t;
        
      } else {
        
        int old_first_parent = data_vertex_firstparent_tetrahedron[v];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_tetrahedra );
        assert( data_tetrahedron_nextparents_of_vertices[ t ][ vi ] == nullindex );
        
        data_vertex_firstparent_tetrahedron[v] = t;
        data_tetrahedron_nextparents_of_vertices[ t ][ vi ] = old_first_parent;
        
      }
      
      assert( data_vertex_firstparent_tetrahedron[v] != nullindex );
      assert( 0 <= data_vertex_firstparent_tetrahedron[v] && data_vertex_firstparent_tetrahedron[v] < counter_tetrahedra );  
    }
    
    
    /* 6. For each vertex, set the first parent face and the neighboring parent faces */
    
    for( int f =  0; f  < counter_faces; f++  )
    for( int vi = 0; vi <             3; vi++ )
    {
      int v = data_face_vertices[f][vi];
      
      assert( 0 <= v && v < counter_vertices );
      
      if( data_vertex_firstparent_face[v] == nullindex ) {
        
        data_vertex_firstparent_face[v] = f;
        
      } else {
        
        int old_first_parent = data_vertex_firstparent_face[v];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_faces );
        assert( data_face_nextparents_of_vertices[ f ][ vi ] == nullindex );
        
        data_vertex_firstparent_face[v] = f;
        data_face_nextparents_of_vertices[ f ][ vi ] = old_first_parent;
        
      }
      
      assert( data_vertex_firstparent_face[v] != nullindex );
      assert( 0 <= data_vertex_firstparent_face[v] && data_vertex_firstparent_face[v] < counter_faces );  
    }
    
    
    /* 7. For each vertex, set the first parent edge and the neighboring parent edges */
    
    for( int e =  0; e  < counter_edges; e++  )
    for( int vi = 0; vi <             2; vi++ )
    {
      int v = data_edge_vertices[e][vi];
      
      assert( 0 <= v && v < counter_vertices );
      
      if( data_vertex_firstparent_edge[v] == nullindex ) {
        
        data_vertex_firstparent_edge[v] = e;
        
      } else {
        
        int old_first_parent = data_vertex_firstparent_edge[v];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_edges );
        assert( data_edge_nextparents_of_vertices[ e ][ vi ] == nullindex );
        
        data_vertex_firstparent_edge[v] = e;
        data_edge_nextparents_of_vertices[ e ][ vi ] = old_first_parent;
        
      }
      
      assert( data_vertex_firstparent_edge[v] != nullindex );
      assert( 0 <= data_vertex_firstparent_edge[v] && data_vertex_firstparent_edge[v] < counter_edges );  
    }
    
    /* 8. For each edge, set the first parent tetrahedron and the neighboring parent tetrahedra */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    for( int e = 0; e <      counter_edges; e++ )
    {
      int voe0 = data_edge_vertices[e][0];
      int voe1 = data_edge_vertices[e][1];
      
      int vot0 = data_tetrahedron_vertices[t][0];
      int vot1 = data_tetrahedron_vertices[t][1];
      int vot2 = data_tetrahedron_vertices[t][2];
      int vot3 = data_tetrahedron_vertices[t][3];
      
      int eot = nullindex;
      
      if( voe0 == vot0 && voe1 == vot1 ) eot = 0;
      if( voe0 == vot0 && voe1 == vot2 ) eot = 1;
      if( voe0 == vot0 && voe1 == vot3 ) eot = 2;
      if( voe0 == vot1 && voe1 == vot2 ) eot = 3;
      if( voe0 == vot1 && voe1 == vot3 ) eot = 4;
      if( voe0 == vot2 && voe1 == vot3 ) eot = 5;
      
      if( eot == nullindex ) continue;
      
      data_tetrahedron_edges[t][eot] = e;
      
      if( data_edge_firstparent_tetrahedron[e] == nullindex ) {
        
        data_edge_firstparent_tetrahedron[e] = t;
        
      } else {
        
        int old_first_parent = data_edge_firstparent_tetrahedron[e];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_tetrahedra );
        
        data_edge_firstparent_tetrahedron[e] = t;
        
        assert( data_tetrahedron_nextparents_of_edges[ t ][ eot ] == nullindex );
        data_tetrahedron_nextparents_of_edges[ t ][ eot ] = old_first_parent;
        
      }
      
      assert( data_edge_firstparent_tetrahedron[e] != nullindex );
      assert( 0 <= data_edge_firstparent_tetrahedron[e] && data_edge_firstparent_tetrahedron[e] < counter_tetrahedra );  
    }
    
    
    /* 9. For each edge, set the first parent face and the neighboring parent faces */
    
    for( int f = 0; f < counter_faces; f++ )
    for( int e = 0; e < counter_edges; e++ )
    {
      int voe0 = data_edge_vertices[e][0];
      int voe1 = data_edge_vertices[e][1];
      
      int vof0 = data_face_vertices[f][0];
      int vof1 = data_face_vertices[f][1];
      int vof2 = data_face_vertices[f][2];
      
      int eof = nullindex;
      
      if( voe0 == vof0 && voe1 == vof1 ) eof = 0;
      if( voe0 == vof0 && voe1 == vof2 ) eof = 1;
      if( voe0 == vof1 && voe1 == vof2 ) eof = 2;
      
      if( eof == nullindex ) continue;
      
      data_face_edges[f][eof] = e;
      
      if( data_edge_firstparent_face[e] == nullindex ) {
        
        data_edge_firstparent_face[e] = f;
        
      } else {
        
        int old_first_parent = data_edge_firstparent_face[e];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_faces );
        
        data_edge_firstparent_face[e] = f;
        
        assert( data_face_nextparents_of_edges[ f ][ eof ] == nullindex );
        data_face_nextparents_of_edges[ f ][ eof ] = old_first_parent;
        
      }
      
      assert( data_edge_firstparent_face[e] != nullindex );
      assert( 0 <= data_edge_firstparent_face[e] && data_edge_firstparent_face[e] < counter_faces );  
    }
    
    
    /* 10. For each face, set the first parent tetrahedron and the neighboring parent tetrahedra */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    for( int f = 0; f <      counter_faces; f++ )
    {
      
      int vof0 = data_face_vertices[f][0];
      int vof1 = data_face_vertices[f][1];
      int vof2 = data_face_vertices[f][2];
      
      int vot0 = data_tetrahedron_vertices[t][0];
      int vot1 = data_tetrahedron_vertices[t][1];
      int vot2 = data_tetrahedron_vertices[t][2];
      int vot3 = data_tetrahedron_vertices[t][3];
      
      
      int fot = nullindex;
      
      if( vof0 == vot0 && vof1 == vot1 && vof2 == vot2 ) fot = 0;
      if( vof0 == vot0 && vof1 == vot1 && vof2 == vot3 ) fot = 1;
      if( vof0 == vot0 && vof1 == vot2 && vof2 == vot3 ) fot = 2;
      if( vof0 == vot1 && vof1 == vot2 && vof2 == vot3 ) fot = 3;
      
      if( fot == nullindex ) continue;
      
      data_tetrahedron_faces[t][fot] = f;
      
      if( data_face_firstparent_tetrahedron[f] == nullindex ) {
        
        data_face_firstparent_tetrahedron[f] = t;
        
      } else {
        
        int old_first_parent = data_face_firstparent_tetrahedron[f];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_tetrahedra );
        
        data_face_firstparent_tetrahedron[f] = t;
        
        assert( data_tetrahedron_nextparents_of_faces[ t ][ fot ] == nullindex );
        data_tetrahedron_nextparents_of_faces[ t ][ fot ] = old_first_parent;
        
      }
      
      assert( data_face_firstparent_tetrahedron[f] != nullindex );
      assert( 0 <= data_face_firstparent_tetrahedron[f] && data_face_firstparent_tetrahedron[f] < counter_tetrahedra );  
    }
    
    
    
    
    check();
}


MeshSimplicial3D::MeshSimplicial3D( 
    int outerdim,
    const Coordinates& coords,
    
    const std::vector<std::array<int,4>> tetrahedron_faces,
    const std::vector<int              > face_firstparent_tetrahedron,
    const std::vector<std::array<int,4>> tetrahedron_nextparents_of_faces,
    
    const std::vector<std::array<int,6>> tetrahedron_edges,
    const std::vector<int              > edge_firstparent_tetrahedron,
    const std::vector<std::array<int,6>> tetrahedron_nextparents_of_edges,
    
    const std::vector<std::array<int,4>> tetrahedron_vertices,
    const std::vector<int              > vertex_firstparent_tetrahedron,
    const std::vector<std::array<int,4>> tetrahedron_nextparents_of_vertices,
    
    const std::vector<std::array<int,3>> face_edges,
    const std::vector<int              > edge_firstparent_face,
    const std::vector<std::array<int,3>> face_nextparents_of_edges,
    
    const std::vector<std::array<int,3>> face_vertices,
    const std::vector<int              > vertex_firstparent_face,
    const std::vector<std::array<int,3>> face_nextparents_of_vertices,
    
    const std::vector<std::array<int,2>> edge_vertices,
    const std::vector<int              > vertex_firstparent_edge,
    const std::vector<std::array<int,2>> edge_nextparents_of_vertices    
)
:
    Mesh( 3, outerdim ),
    
    counter_tetrahedra( tetrahedron_vertices.size() ),
    counter_faces( face_vertices.size() ),
    counter_edges( edge_vertices.size() ),
    counter_vertices( vertex_firstparent_tetrahedron.size() ),
    
    data_tetrahedron_faces( tetrahedron_faces ),
    data_face_firstparent_tetrahedron( face_firstparent_tetrahedron ),
    data_tetrahedron_nextparents_of_faces( tetrahedron_nextparents_of_faces ),
    
    data_tetrahedron_edges( tetrahedron_edges ),
    data_edge_firstparent_tetrahedron( edge_firstparent_tetrahedron ),
    data_tetrahedron_nextparents_of_edges( tetrahedron_nextparents_of_edges ),
    
    data_tetrahedron_vertices( tetrahedron_vertices ),
    data_vertex_firstparent_tetrahedron( vertex_firstparent_tetrahedron ),
    data_tetrahedron_nextparents_of_vertices( tetrahedron_nextparents_of_vertices ),
    
    data_face_edges( face_edges ),
    data_edge_firstparent_face( edge_firstparent_face ),
    data_face_nextparents_of_edges( face_nextparents_of_edges ),
    
    data_face_vertices( face_vertices ),
    data_vertex_firstparent_face( vertex_firstparent_face ),
    data_face_nextparents_of_vertices( face_nextparents_of_vertices ),
    
    data_edge_vertices( edge_vertices ),
    data_vertex_firstparent_edge( vertex_firstparent_edge ),
    data_edge_nextparents_of_vertices( edge_nextparents_of_vertices )
{
    
    getcoordinates() = coords;
    
    check();
}


MeshSimplicial3D::~MeshSimplicial3D()
{
    
}

bool MeshSimplicial3D::operator== ( const MeshSimplicial3D& mesh ) const 
{
  return counter_tetrahedra == mesh.counter_tetrahedra
         &&
         counter_faces == mesh.counter_faces
         &&
         counter_edges == mesh.counter_edges
         &&
         counter_vertices == mesh.counter_vertices
         &&
         data_tetrahedron_faces == mesh.data_tetrahedron_faces
         &&
         data_face_firstparent_tetrahedron == mesh.data_face_firstparent_tetrahedron
         &&
         data_tetrahedron_nextparents_of_faces == mesh.data_tetrahedron_nextparents_of_faces
         &&
         data_tetrahedron_edges == mesh.data_tetrahedron_edges
         &&
         data_edge_firstparent_tetrahedron == mesh.data_edge_firstparent_tetrahedron
         &&
         data_tetrahedron_nextparents_of_edges == mesh.data_tetrahedron_nextparents_of_edges
         &&
         data_tetrahedron_vertices == mesh.data_tetrahedron_vertices 
         &&
         data_vertex_firstparent_tetrahedron == mesh.data_vertex_firstparent_tetrahedron
         &&
         data_tetrahedron_nextparents_of_vertices == mesh.data_tetrahedron_nextparents_of_vertices
         &&
         data_face_edges == mesh.data_face_edges
         &&
         data_edge_firstparent_face == mesh.data_edge_firstparent_face
         &&
         data_face_nextparents_of_edges == mesh.data_face_nextparents_of_edges
         &&
         data_face_vertices == mesh.data_face_vertices 
         &&
         data_vertex_firstparent_face == mesh.data_vertex_firstparent_face
         &&
         data_face_nextparents_of_vertices == mesh.data_face_nextparents_of_vertices
         &&
         data_edge_vertices == mesh.data_edge_vertices 
         &&
         data_vertex_firstparent_edge == mesh.data_vertex_firstparent_edge
         &&
         data_edge_nextparents_of_vertices == mesh.data_edge_nextparents_of_vertices
         &&
         getinnerdimension() == mesh.getinnerdimension()
         &&
         getouterdimension() == mesh.getouterdimension()
         &&
         getcoordinates() == mesh.getcoordinates()
         &&
         true;
}

bool MeshSimplicial3D::operator!= ( const MeshSimplicial3D& mesh ) const 
{
  return ! ( *this == mesh );
}



















void MeshSimplicial3D::check() const
{
    
    /****************************/
    /* 1. Check the array sizes */ // OK
    /****************************/
    
    assert( counter_tetrahedra == data_tetrahedron_faces.size() );
    assert( counter_tetrahedra == data_tetrahedron_nextparents_of_faces.size() );
    assert( counter_tetrahedra == data_tetrahedron_edges.size() );
    assert( counter_tetrahedra == data_tetrahedron_nextparents_of_edges.size() );
    assert( counter_tetrahedra == data_tetrahedron_vertices.size() );
    assert( counter_tetrahedra == data_tetrahedron_nextparents_of_vertices.size() );
    
    assert( counter_faces == data_face_edges.size() );
    assert( counter_faces == data_face_nextparents_of_edges.size() );
    assert( counter_faces == data_face_vertices.size() );
    assert( counter_faces == data_face_nextparents_of_vertices.size() );
    
    assert( counter_faces == data_face_firstparent_tetrahedron.size() );
    
    assert( counter_edges == data_edge_vertices.size() );
    assert( counter_edges == data_edge_nextparents_of_vertices.size() );
    
    assert( counter_edges == data_edge_firstparent_face.size() );
    assert( counter_edges == data_edge_firstparent_tetrahedron.size() );
    
    assert( counter_vertices == data_vertex_firstparent_edge.size() );
    assert( counter_vertices == data_vertex_firstparent_face.size() );
    assert( counter_vertices == data_vertex_firstparent_tetrahedron.size() );
    
    assert( count_vertices() == getcoordinates().getnumber() );
    
    
    
    
    /********************************************************************************/
    /* 2. check that the internal data of each simplex make sense on each dimension */
    /********************************************************************************/
    
    /* VERTEX EDGE STUFF */
    
    /* each edge: each vertex is a valid index          */
    /* each edge: each vertex is unique                 */
    
    for( int e = 0; e < counter_edges; e++ )
    {

        assert( data_edge_vertices[e][0] != nullindex );
        assert( data_edge_vertices[e][1] != nullindex );
        assert( 0 <= data_edge_vertices[e][0] && data_edge_vertices[e][0] < counter_vertices );
        assert( 0 <= data_edge_vertices[e][1] && data_edge_vertices[e][0] < counter_vertices );
        assert( data_edge_vertices[e][0] != data_edge_vertices[e][1] );
        
    }
    
    /* VERTEX FACE STUFF */
    
    /* each face: each vertex is a valid index          */
    /* each face: each vertex is unique                 */
    
    for( int f = 0; f < counter_faces; f++ )
    {
        
        assert( data_face_vertices[f][0] != nullindex );
        assert( data_face_vertices[f][1] != nullindex );
        assert( data_face_vertices[f][2] != nullindex );
        
        assert( 0 <= data_face_vertices[f][0] && data_face_vertices[f][0] < counter_vertices );
        assert( 0 <= data_face_vertices[f][1] && data_face_vertices[f][1] < counter_vertices );
        assert( 0 <= data_face_vertices[f][2] && data_face_vertices[f][2] < counter_vertices );
        
        assert( data_face_vertices[f][0] != data_face_vertices[f][1] );
        assert( data_face_vertices[f][0] != data_face_vertices[f][2] );
        assert( data_face_vertices[f][1] != data_face_vertices[f][2] );
       
    }
    
    /* EDGE FACE STUFF */
    
    /* each face: each edge is a valid index            */
    /* each face: each edge is unique                   */
    
    for( int f = 0; f < counter_faces; f++ )
    {
        
        assert( data_face_edges[f][0] != nullindex );
        assert( data_face_edges[f][1] != nullindex );
        assert( data_face_edges[f][2] != nullindex );
        
        assert( 0 <= data_face_edges[f][0] && data_face_edges[f][0] < counter_edges );
        assert( 0 <= data_face_edges[f][1] && data_face_edges[f][1] < counter_edges );
        assert( 0 <= data_face_edges[f][2] && data_face_edges[f][2] < counter_edges );
        
        assert( data_face_edges[f][0] != data_face_edges[f][1] );
        assert( data_face_edges[f][0] != data_face_edges[f][2] );
        assert( data_face_edges[f][1] != data_face_edges[f][2] );
        
    }
    
    /* VERTEX TET STUFF */
    
    /* each tet: each vertex is a valid index */
    /* each tet: each vertex is unique        */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        
        assert( data_tetrahedron_vertices[t][0] != nullindex );
        assert( data_tetrahedron_vertices[t][1] != nullindex );
        assert( data_tetrahedron_vertices[t][2] != nullindex );
        assert( data_tetrahedron_vertices[t][3] != nullindex );
        
        assert( 0 <= data_tetrahedron_vertices[t][0] && data_tetrahedron_vertices[t][0] < counter_vertices );
        assert( 0 <= data_tetrahedron_vertices[t][1] && data_tetrahedron_vertices[t][1] < counter_vertices );
        assert( 0 <= data_tetrahedron_vertices[t][2] && data_tetrahedron_vertices[t][2] < counter_vertices );
        assert( 0 <= data_tetrahedron_vertices[t][3] && data_tetrahedron_vertices[t][3] < counter_vertices );
        
        assert( data_tetrahedron_vertices[t][0] != data_tetrahedron_vertices[t][1] );
        assert( data_tetrahedron_vertices[t][0] != data_tetrahedron_vertices[t][2] );
        assert( data_tetrahedron_vertices[t][0] != data_tetrahedron_vertices[t][3] );
        assert( data_tetrahedron_vertices[t][1] != data_tetrahedron_vertices[t][2] );
        assert( data_tetrahedron_vertices[t][1] != data_tetrahedron_vertices[t][3] );
        assert( data_tetrahedron_vertices[t][2] != data_tetrahedron_vertices[t][3] );
        
    }
    
    /* EDGE TET STUFF */
    
    /* each tet: each edge is a valid index */
    /* each tet: each edge is unique        */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        
        assert( data_tetrahedron_edges[t][0] != nullindex );
        assert( data_tetrahedron_edges[t][1] != nullindex );
        assert( data_tetrahedron_edges[t][2] != nullindex );
        assert( data_tetrahedron_edges[t][3] != nullindex );
        assert( data_tetrahedron_edges[t][4] != nullindex );
        assert( data_tetrahedron_edges[t][5] != nullindex );
        
        assert( 0 <= data_tetrahedron_edges[t][0] && data_tetrahedron_edges[t][0] < counter_edges );
        assert( 0 <= data_tetrahedron_edges[t][1] && data_tetrahedron_edges[t][1] < counter_edges );
        assert( 0 <= data_tetrahedron_edges[t][2] && data_tetrahedron_edges[t][2] < counter_edges );
        assert( 0 <= data_tetrahedron_edges[t][3] && data_tetrahedron_edges[t][3] < counter_edges );
        assert( 0 <= data_tetrahedron_edges[t][4] && data_tetrahedron_edges[t][4] < counter_edges );
        assert( 0 <= data_tetrahedron_edges[t][5] && data_tetrahedron_edges[t][5] < counter_edges );
        
        assert( data_tetrahedron_edges[t][0] != data_tetrahedron_edges[t][1] );
        assert( data_tetrahedron_edges[t][0] != data_tetrahedron_edges[t][2] );
        assert( data_tetrahedron_edges[t][0] != data_tetrahedron_edges[t][3] );
        assert( data_tetrahedron_edges[t][0] != data_tetrahedron_edges[t][4] );
        assert( data_tetrahedron_edges[t][0] != data_tetrahedron_edges[t][5] );
        assert( data_tetrahedron_edges[t][1] != data_tetrahedron_edges[t][2] );
        assert( data_tetrahedron_edges[t][1] != data_tetrahedron_edges[t][3] );
        assert( data_tetrahedron_edges[t][1] != data_tetrahedron_edges[t][4] );
        assert( data_tetrahedron_edges[t][1] != data_tetrahedron_edges[t][5] );
        assert( data_tetrahedron_edges[t][2] != data_tetrahedron_edges[t][3] );
        assert( data_tetrahedron_edges[t][2] != data_tetrahedron_edges[t][4] );
        assert( data_tetrahedron_edges[t][2] != data_tetrahedron_edges[t][5] );
        assert( data_tetrahedron_edges[t][3] != data_tetrahedron_edges[t][4] );
        assert( data_tetrahedron_edges[t][3] != data_tetrahedron_edges[t][5] );
        assert( data_tetrahedron_edges[t][4] != data_tetrahedron_edges[t][5] );
        
    }
    
    /* FACE TET STUFF */
    
    /* each tet: each face is a valid index */
    /* each tet: each face is unique        */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        
        assert( data_tetrahedron_faces[t][0] != nullindex );
        assert( data_tetrahedron_faces[t][1] != nullindex );
        assert( data_tetrahedron_faces[t][2] != nullindex );
        assert( data_tetrahedron_faces[t][3] != nullindex );
        
        assert( 0 <= data_tetrahedron_faces[t][0] && data_tetrahedron_faces[t][0] < counter_faces );
        assert( 0 <= data_tetrahedron_faces[t][1] && data_tetrahedron_faces[t][1] < counter_faces );
        assert( 0 <= data_tetrahedron_faces[t][2] && data_tetrahedron_faces[t][2] < counter_faces );
        assert( 0 <= data_tetrahedron_faces[t][3] && data_tetrahedron_faces[t][3] < counter_faces );
        
        assert( data_tetrahedron_faces[t][0] != data_tetrahedron_faces[t][1] );
        assert( data_tetrahedron_faces[t][0] != data_tetrahedron_faces[t][2] );
        assert( data_tetrahedron_faces[t][0] != data_tetrahedron_faces[t][3] );
        assert( data_tetrahedron_faces[t][1] != data_tetrahedron_faces[t][2] );
        assert( data_tetrahedron_faces[t][1] != data_tetrahedron_faces[t][3] );
        assert( data_tetrahedron_faces[t][2] != data_tetrahedron_faces[t][3] );
        
    }
    
    
    
    
    /****************************************/
    /* 2. check the uniqueness of simplices */
    /****************************************/
    
    /* 
     * check that all edges are unique, even up to permutation
     */
    
    for( int e1 = 0; e1 < counter_edges; e1++ )
    for( int e2 = 0; e2 < counter_edges; e2++ )
    {
        if( e1 == e2 ) continue;
        
        assert( data_edge_vertices[e1][0] != data_edge_vertices[e2][0] || data_edge_vertices[e1][1] != data_edge_vertices[e2][1] );
        assert( data_edge_vertices[e1][0] != data_edge_vertices[e2][1] || data_edge_vertices[e1][1] != data_edge_vertices[e2][0] );
    }
    
    /* 
     * check that all faces are unique, even up to permutation
     * check both the vertices and the edges listed for each face
     */
    
    for( int f1 = 0; f1 < counter_faces; f1++ )
    for( int f2 = 0; f2 < counter_faces; f2++ )
    {
        if( f1 == f2 ) continue;
        
        assert( data_face_vertices[f1][0] != data_face_vertices[f2][0] || data_face_vertices[f1][1] != data_face_vertices[f2][1] || data_face_vertices[f1][2] != data_face_vertices[f2][2] );
        assert( data_face_vertices[f1][0] != data_face_vertices[f2][0] || data_face_vertices[f1][1] != data_face_vertices[f2][2] || data_face_vertices[f1][2] != data_face_vertices[f2][1] );
        assert( data_face_vertices[f1][0] != data_face_vertices[f2][1] || data_face_vertices[f1][1] != data_face_vertices[f2][0] || data_face_vertices[f1][2] != data_face_vertices[f2][2] );
        assert( data_face_vertices[f1][0] != data_face_vertices[f2][1] || data_face_vertices[f1][1] != data_face_vertices[f2][2] || data_face_vertices[f1][2] != data_face_vertices[f2][0] );
        assert( data_face_vertices[f1][0] != data_face_vertices[f2][2] || data_face_vertices[f1][1] != data_face_vertices[f2][0] || data_face_vertices[f1][2] != data_face_vertices[f2][1] );
        assert( data_face_vertices[f1][0] != data_face_vertices[f2][2] || data_face_vertices[f1][1] != data_face_vertices[f2][1] || data_face_vertices[f1][2] != data_face_vertices[f2][0] );
        
        assert( data_face_edges[f1][0] != data_face_edges[f2][0] || data_face_edges[f1][1] != data_face_edges[f2][1] || data_face_edges[f1][2] != data_face_edges[f2][2] );
        assert( data_face_edges[f1][0] != data_face_edges[f2][0] || data_face_edges[f1][1] != data_face_edges[f2][2] || data_face_edges[f1][2] != data_face_edges[f2][1] );
        assert( data_face_edges[f1][0] != data_face_edges[f2][1] || data_face_edges[f1][1] != data_face_edges[f2][0] || data_face_edges[f1][2] != data_face_edges[f2][2] );
        assert( data_face_edges[f1][0] != data_face_edges[f2][1] || data_face_edges[f1][1] != data_face_edges[f2][2] || data_face_edges[f1][2] != data_face_edges[f2][0] );
        assert( data_face_edges[f1][0] != data_face_edges[f2][2] || data_face_edges[f1][1] != data_face_edges[f2][0] || data_face_edges[f1][2] != data_face_edges[f2][1] );
        assert( data_face_edges[f1][0] != data_face_edges[f2][2] || data_face_edges[f1][1] != data_face_edges[f2][1] || data_face_edges[f1][2] != data_face_edges[f2][0] );
        
    }
    
    /* 
     * check that all tetrahedra are unique, even up to permutation
     * check both the vertices, the edges, and the faces listed for each tetrahedron
     */
    
    const auto permutation_of_six_objects = generatePermutations( IndexRange( 0, 5 ) );
      
    for( int t1 = 0; t1 < counter_tetrahedra; t1++ )
    for( int t2 = 0; t2 < counter_tetrahedra; t2++ )
    {
        if( t1 == t2 ) continue;
        
        /* vertices */
        
        const auto& refv = data_tetrahedron_vertices;
        assert( refv[t1][0] != refv[t2][0] || refv[t1][1] != refv[t2][1] || refv[t1][2] != refv[t2][2] || refv[t1][3] != refv[t2][3] );
        assert( refv[t1][0] != refv[t2][0] || refv[t1][1] != refv[t2][1] || refv[t1][2] != refv[t2][3] || refv[t1][3] != refv[t2][2] );
        assert( refv[t1][0] != refv[t2][0] || refv[t1][1] != refv[t2][2] || refv[t1][2] != refv[t2][1] || refv[t1][3] != refv[t2][3] );
        assert( refv[t1][0] != refv[t2][0] || refv[t1][1] != refv[t2][2] || refv[t1][2] != refv[t2][3] || refv[t1][3] != refv[t2][1] );
        assert( refv[t1][0] != refv[t2][0] || refv[t1][1] != refv[t2][3] || refv[t1][2] != refv[t2][1] || refv[t1][3] != refv[t2][2] );
        assert( refv[t1][0] != refv[t2][0] || refv[t1][1] != refv[t2][3] || refv[t1][2] != refv[t2][2] || refv[t1][3] != refv[t2][1] );
        
        assert( refv[t1][0] != refv[t2][1] || refv[t1][1] != refv[t2][0] || refv[t1][2] != refv[t2][2] || refv[t1][3] != refv[t2][3] );
        assert( refv[t1][0] != refv[t2][1] || refv[t1][1] != refv[t2][0] || refv[t1][2] != refv[t2][1] || refv[t1][3] != refv[t2][2] );
        assert( refv[t1][0] != refv[t2][1] || refv[t1][1] != refv[t2][2] || refv[t1][2] != refv[t2][0] || refv[t1][3] != refv[t2][3] );
        assert( refv[t1][0] != refv[t2][1] || refv[t1][1] != refv[t2][2] || refv[t1][2] != refv[t2][3] || refv[t1][3] != refv[t2][0] );
        assert( refv[t1][0] != refv[t2][1] || refv[t1][1] != refv[t2][3] || refv[t1][2] != refv[t2][0] || refv[t1][3] != refv[t2][2] );
        assert( refv[t1][0] != refv[t2][1] || refv[t1][1] != refv[t2][3] || refv[t1][2] != refv[t2][2] || refv[t1][3] != refv[t2][0] );
        
        assert( refv[t1][0] != refv[t2][2] || refv[t1][1] != refv[t2][0] || refv[t1][2] != refv[t2][1] || refv[t1][3] != refv[t2][3] );
        assert( refv[t1][0] != refv[t2][2] || refv[t1][1] != refv[t2][0] || refv[t1][2] != refv[t2][3] || refv[t1][3] != refv[t2][1] );
        assert( refv[t1][0] != refv[t2][2] || refv[t1][1] != refv[t2][1] || refv[t1][2] != refv[t2][0] || refv[t1][3] != refv[t2][3] );
        assert( refv[t1][0] != refv[t2][2] || refv[t1][1] != refv[t2][1] || refv[t1][2] != refv[t2][3] || refv[t1][3] != refv[t2][0] );
        assert( refv[t1][0] != refv[t2][2] || refv[t1][1] != refv[t2][3] || refv[t1][2] != refv[t2][0] || refv[t1][3] != refv[t2][1] );
        assert( refv[t1][0] != refv[t2][2] || refv[t1][1] != refv[t2][3] || refv[t1][2] != refv[t2][1] || refv[t1][3] != refv[t2][0] );
        
        assert( refv[t1][0] != refv[t2][3] || refv[t1][1] != refv[t2][0] || refv[t1][2] != refv[t2][1] || refv[t1][3] != refv[t2][2] );
        assert( refv[t1][0] != refv[t2][3] || refv[t1][1] != refv[t2][0] || refv[t1][2] != refv[t2][2] || refv[t1][3] != refv[t2][1] );
        assert( refv[t1][0] != refv[t2][3] || refv[t1][1] != refv[t2][1] || refv[t1][2] != refv[t2][0] || refv[t1][3] != refv[t2][2] );
        assert( refv[t1][0] != refv[t2][3] || refv[t1][1] != refv[t2][1] || refv[t1][2] != refv[t2][2] || refv[t1][3] != refv[t2][0] );
        assert( refv[t1][0] != refv[t2][3] || refv[t1][1] != refv[t2][2] || refv[t1][2] != refv[t2][0] || refv[t1][3] != refv[t2][1] );
        assert( refv[t1][0] != refv[t2][3] || refv[t1][1] != refv[t2][2] || refv[t1][2] != refv[t2][1] || refv[t1][3] != refv[t2][0] );
        
        
        /* edges */
        
        const auto& refe = data_tetrahedron_edges;
        
        const auto& perms = permutation_of_six_objects;
      
        for( auto perm : perms )
          assert( refe[t1][0] != refe[t2][ perm[0] ] || refe[t1][1] != refe[t2][ perm[1] ] || refe[t1][2] != refe[t2][ perm[2] ]
                  ||
                  refe[t1][3] != refe[t2][ perm[3] ] || refe[t1][4] != refe[t2][ perm[4] ] || refe[t1][5] != refe[t2][ perm[5] ] );
        
        
        /* faces */
        
        const auto& reff = data_tetrahedron_faces;
        assert( reff[t1][0] != reff[t2][0] || reff[t1][1] != reff[t2][1] || reff[t1][2] != reff[t2][2] || reff[t1][3] != reff[t2][3] );
        assert( reff[t1][0] != reff[t2][0] || reff[t1][1] != reff[t2][1] || reff[t1][2] != reff[t2][3] || reff[t1][3] != reff[t2][2] );
        assert( reff[t1][0] != reff[t2][0] || reff[t1][1] != reff[t2][2] || reff[t1][2] != reff[t2][1] || reff[t1][3] != reff[t2][3] );
        assert( reff[t1][0] != reff[t2][0] || reff[t1][1] != reff[t2][2] || reff[t1][2] != reff[t2][3] || reff[t1][3] != reff[t2][1] );
        assert( reff[t1][0] != reff[t2][0] || reff[t1][1] != reff[t2][3] || reff[t1][2] != reff[t2][1] || reff[t1][3] != reff[t2][2] );
        assert( reff[t1][0] != reff[t2][0] || reff[t1][1] != reff[t2][3] || reff[t1][2] != reff[t2][2] || reff[t1][3] != reff[t2][1] );
        
        assert( reff[t1][0] != reff[t2][1] || reff[t1][1] != reff[t2][0] || reff[t1][2] != reff[t2][2] || reff[t1][3] != reff[t2][3] );
        assert( reff[t1][0] != reff[t2][1] || reff[t1][1] != reff[t2][0] || reff[t1][2] != reff[t2][1] || reff[t1][3] != reff[t2][2] );
        assert( reff[t1][0] != reff[t2][1] || reff[t1][1] != reff[t2][2] || reff[t1][2] != reff[t2][0] || reff[t1][3] != reff[t2][3] );
        assert( reff[t1][0] != reff[t2][1] || reff[t1][1] != reff[t2][2] || reff[t1][2] != reff[t2][3] || reff[t1][3] != reff[t2][0] );
        assert( reff[t1][0] != reff[t2][1] || reff[t1][1] != reff[t2][3] || reff[t1][2] != reff[t2][0] || reff[t1][3] != reff[t2][2] );
        assert( reff[t1][0] != reff[t2][1] || reff[t1][1] != reff[t2][3] || reff[t1][2] != reff[t2][2] || reff[t1][3] != reff[t2][0] );
        
        assert( reff[t1][0] != reff[t2][2] || reff[t1][1] != reff[t2][0] || reff[t1][2] != reff[t2][1] || reff[t1][3] != reff[t2][3] );
        assert( reff[t1][0] != reff[t2][2] || reff[t1][1] != reff[t2][0] || reff[t1][2] != reff[t2][3] || reff[t1][3] != reff[t2][1] );
        assert( reff[t1][0] != reff[t2][2] || reff[t1][1] != reff[t2][1] || reff[t1][2] != reff[t2][0] || reff[t1][3] != reff[t2][3] );
        assert( reff[t1][0] != reff[t2][2] || reff[t1][1] != reff[t2][1] || reff[t1][2] != reff[t2][3] || reff[t1][3] != reff[t2][0] );
        assert( reff[t1][0] != reff[t2][2] || reff[t1][1] != reff[t2][3] || reff[t1][2] != reff[t2][0] || reff[t1][3] != reff[t2][1] );
        assert( reff[t1][0] != reff[t2][2] || reff[t1][1] != reff[t2][3] || reff[t1][2] != reff[t2][1] || reff[t1][3] != reff[t2][0] );
        
        assert( reff[t1][0] != reff[t2][3] || reff[t1][1] != reff[t2][0] || reff[t1][2] != reff[t2][1] || reff[t1][3] != reff[t2][2] );
        assert( reff[t1][0] != reff[t2][3] || reff[t1][1] != reff[t2][0] || reff[t1][2] != reff[t2][2] || reff[t1][3] != reff[t2][1] );
        assert( reff[t1][0] != reff[t2][3] || reff[t1][1] != reff[t2][1] || reff[t1][2] != reff[t2][0] || reff[t1][3] != reff[t2][2] );
        assert( reff[t1][0] != reff[t2][3] || reff[t1][1] != reff[t2][1] || reff[t1][2] != reff[t2][2] || reff[t1][3] != reff[t2][0] );
        assert( reff[t1][0] != reff[t2][3] || reff[t1][1] != reff[t2][2] || reff[t1][2] != reff[t2][0] || reff[t1][3] != reff[t2][1] );
        assert( reff[t1][0] != reff[t2][3] || reff[t1][1] != reff[t2][2] || reff[t1][2] != reff[t2][1] || reff[t1][3] != reff[t2][0] );
        
    }
    
    
    
    
    
    
    
    /********************************************************/
    /* 3. check the data of each simplex accross dimensions */
    /********************************************************/
    
    
    /*
     * each face: each edge is listed correctly
     */
    
    for( int f = 0; f < counter_faces; f++ )
    {
        assert( data_edge_vertices[ data_face_edges[f][0] ][0] == data_face_vertices[f][0] );
        assert( data_edge_vertices[ data_face_edges[f][0] ][1] == data_face_vertices[f][1] );
        
        assert( data_edge_vertices[ data_face_edges[f][1] ][0] == data_face_vertices[f][0] );
        assert( data_edge_vertices[ data_face_edges[f][1] ][1] == data_face_vertices[f][2] );
        
        assert( data_edge_vertices[ data_face_edges[f][2] ][0] == data_face_vertices[f][1] );
        assert( data_edge_vertices[ data_face_edges[f][2] ][1] == data_face_vertices[f][2] );
    }
    
    
    
    /*
     * each tetrahedron: each edge is listed correctly
     */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        assert( data_edge_vertices[ data_tetrahedron_edges[t][0] ][0] == data_tetrahedron_vertices[t][0] );
        assert( data_edge_vertices[ data_tetrahedron_edges[t][0] ][1] == data_tetrahedron_vertices[t][1] );
        
        assert( data_edge_vertices[ data_tetrahedron_edges[t][1] ][0] == data_tetrahedron_vertices[t][0] );
        assert( data_edge_vertices[ data_tetrahedron_edges[t][1] ][1] == data_tetrahedron_vertices[t][2] );
        
        assert( data_edge_vertices[ data_tetrahedron_edges[t][2] ][0] == data_tetrahedron_vertices[t][0] );
        assert( data_edge_vertices[ data_tetrahedron_edges[t][2] ][1] == data_tetrahedron_vertices[t][3] );
        
        assert( data_edge_vertices[ data_tetrahedron_edges[t][3] ][0] == data_tetrahedron_vertices[t][1] );
        assert( data_edge_vertices[ data_tetrahedron_edges[t][3] ][1] == data_tetrahedron_vertices[t][2] );
        
        assert( data_edge_vertices[ data_tetrahedron_edges[t][4] ][0] == data_tetrahedron_vertices[t][1] );
        assert( data_edge_vertices[ data_tetrahedron_edges[t][4] ][1] == data_tetrahedron_vertices[t][3] );
        
        assert( data_edge_vertices[ data_tetrahedron_edges[t][5] ][0] == data_tetrahedron_vertices[t][2] );
        assert( data_edge_vertices[ data_tetrahedron_edges[t][5] ][1] == data_tetrahedron_vertices[t][3] );
    }
    
    
    /*
     * each tetrahedron: each face is listed correctly
     */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        assert( data_face_vertices[ data_tetrahedron_faces[t][0] ][0] == data_tetrahedron_vertices[t][0] );
        assert( data_face_vertices[ data_tetrahedron_faces[t][0] ][1] == data_tetrahedron_vertices[t][1] );
        assert( data_face_vertices[ data_tetrahedron_faces[t][0] ][2] == data_tetrahedron_vertices[t][2] );
        
        assert( data_face_vertices[ data_tetrahedron_faces[t][1] ][0] == data_tetrahedron_vertices[t][0] );
        assert( data_face_vertices[ data_tetrahedron_faces[t][1] ][1] == data_tetrahedron_vertices[t][1] );
        assert( data_face_vertices[ data_tetrahedron_faces[t][1] ][2] == data_tetrahedron_vertices[t][3] );
        
        assert( data_face_vertices[ data_tetrahedron_faces[t][2] ][0] == data_tetrahedron_vertices[t][0] );
        assert( data_face_vertices[ data_tetrahedron_faces[t][2] ][1] == data_tetrahedron_vertices[t][2] );
        assert( data_face_vertices[ data_tetrahedron_faces[t][2] ][2] == data_tetrahedron_vertices[t][3] );
        
        assert( data_face_vertices[ data_tetrahedron_faces[t][3] ][0] == data_tetrahedron_vertices[t][1] );
        assert( data_face_vertices[ data_tetrahedron_faces[t][3] ][1] == data_tetrahedron_vertices[t][2] );
        assert( data_face_vertices[ data_tetrahedron_faces[t][3] ][2] == data_tetrahedron_vertices[t][3] );
    }
    
    
    
    
    /**************************************************************/
    /* 4. Check that the first parents are set and actual parents */
    /**************************************************************/
    
    /* each first parent tetrahedron of an face: first parent is non-null and a valid parent */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        int p = data_face_firstparent_tetrahedron[t];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_tetrahedra );
        
        assert( data_tetrahedron_faces[p][0] == t || data_tetrahedron_faces[p][1] == t || data_tetrahedron_faces[p][2] == t || data_tetrahedron_faces[p][3] == t );
    }
    
    /* each first parent tetrahedron of an edge: first parent is non-null and a valid parent */
    
    for( int e = 0; e < counter_edges; e++ )
    {
        int p = data_edge_firstparent_tetrahedron[e];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_tetrahedra );
        
        assert( data_tetrahedron_edges[p][0] == e || data_tetrahedron_edges[p][1] == e || data_tetrahedron_edges[p][2] == e 
                || 
                data_tetrahedron_edges[p][3] == e || data_tetrahedron_edges[p][4] == e || data_tetrahedron_edges[p][5] == e
              );
    }
    
    /* each first parent tetrahedron of an vertex: first parent is non-null and a valid parent */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        int p = data_vertex_firstparent_tetrahedron[t];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_tetrahedra );
        
        assert( data_tetrahedron_vertices[p][0] == t || data_tetrahedron_vertices[p][1] == t || data_tetrahedron_vertices[p][2] == t || data_tetrahedron_vertices[p][3] == t );
    }
    
    /* each first parent face of an edge: first parent is non-null and a valid parent */
    
    for( int e = 0; e < counter_edges; e++ )
    {
        int p = data_edge_firstparent_face[e];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_faces );
        
        assert( data_face_edges[p][0] == e || data_face_edges[p][1] == e || data_face_edges[p][2] == e );
    }
    
    /* each first parent face of a vertex: first parent is non-null and a valid parent */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
        int p = data_vertex_firstparent_face[v];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_faces );
        
        assert( data_face_vertices[p][0] == v || data_face_vertices[p][1] == v || data_face_vertices[p][2] == v );
    }
    
    /* each first parent edge of a vertex: first parent is non-null and a valid parent */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
        int p = data_vertex_firstparent_edge[v];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_edges );
        
        assert( data_edge_vertices[p][0] == v || data_edge_vertices[p][1] == v );
    }
    
    
    
    
    
    
    
    
    /*****************************************************/
    /* 5. Check that the next parents are actual parents */
    /*****************************************************/
    
    /* VERTEX EDGE STUFF */
    
    /* each edge: the next parents make sense           */ 
    /* each edge: the next parents are actually parents */
    
    for( int e = 0; e < counter_edges; e++ )
    for( int vi = 0; vi < 2; vi++ )
    {
        
        if( data_edge_nextparents_of_vertices[e][vi] != nullindex )
        assert( 0 <= data_edge_nextparents_of_vertices[e][vi] && data_edge_nextparents_of_vertices[e][vi] < counter_edges );
        
        if( data_edge_nextparents_of_vertices[e][vi] != nullindex )
        assert( data_edge_vertices[ data_edge_nextparents_of_vertices[e][vi] ][0] == data_edge_vertices[e][vi] 
                ||
                data_edge_vertices[ data_edge_nextparents_of_vertices[e][vi] ][1] == data_edge_vertices[e][vi] 
                );
        
    }
    
    /* VERTEX FACE STUFF */
    
    /* each face: the next parents are unique           */
    /* each face: the next parents are actually parents */
    
    for( int f = 0; f < counter_faces; f++ )
    for( int vi = 0; vi < 3; vi++ )
    {
        
        if( data_face_nextparents_of_vertices[f][vi] != nullindex )
            assert( 0 <= data_face_nextparents_of_vertices[f][vi] && data_face_nextparents_of_vertices[f][vi] < counter_faces );
    
        if( data_face_nextparents_of_vertices[f][vi] != nullindex )
            assert( data_face_vertices[ data_face_nextparents_of_vertices[f][vi] ][0] == data_face_vertices[f][vi] 
                    ||
                    data_face_vertices[ data_face_nextparents_of_vertices[f][vi] ][1] == data_face_vertices[f][vi] 
                    ||
                    data_face_vertices[ data_face_nextparents_of_vertices[f][vi] ][2] == data_face_vertices[f][vi] );
        
    }
    
    /* EDGE FACE STUFF */
    
    /* each face: the next parents are unique           */
    /* each face: the next parents are actually parents */
    
    for( int f = 0; f < counter_faces; f++ )
    for( int ei = 0; ei < 3; ei++ )
    {
        
        if( data_face_nextparents_of_edges[f][ei] != nullindex )
            assert( 0 <= data_face_nextparents_of_edges[f][ei] && data_face_nextparents_of_edges[f][ei] < counter_faces );
    
        if( data_face_nextparents_of_edges[f][ei] != nullindex )
            assert( data_face_edges[ data_face_nextparents_of_edges[f][ei] ][0] == data_face_edges[f][ei] 
                    ||
                    data_face_edges[ data_face_nextparents_of_edges[f][ei] ][1] == data_face_edges[f][ei] 
                    ||
                    data_face_edges[ data_face_nextparents_of_edges[f][ei] ][2] == data_face_edges[f][ei] );
        
    }
    
    /* VERTEX TET STUFF */
    
    /* each tet: the next parents are unique           */
    /* each tet: the next parents are actually parents */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    for( int vi = 0; vi < 4; vi++ )
    {
        
        if( data_tetrahedron_nextparents_of_vertices[t][vi] != nullindex )
            assert( 0 <= data_tetrahedron_nextparents_of_vertices[t][vi] && data_tetrahedron_nextparents_of_vertices[t][vi] < counter_vertices );
    
        if( data_tetrahedron_nextparents_of_vertices[t][vi] != nullindex )
            assert( data_tetrahedron_vertices[ data_tetrahedron_nextparents_of_vertices[t][vi] ][0] == data_tetrahedron_vertices[t][vi] 
                    ||
                    data_tetrahedron_vertices[ data_tetrahedron_nextparents_of_vertices[t][vi] ][1] == data_tetrahedron_vertices[t][vi] 
                    ||
                    data_tetrahedron_vertices[ data_tetrahedron_nextparents_of_vertices[t][vi] ][2] == data_tetrahedron_vertices[t][vi] 
                    ||
                    data_tetrahedron_vertices[ data_tetrahedron_nextparents_of_vertices[t][vi] ][3] == data_tetrahedron_vertices[t][vi] );
        
    }
    
    /* EDGE TET STUFF */
    
    /* each tet: the next parents are unique           */
    /* each tet: the next parents are actually parents */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    for( int ei = 0; ei < 6; ei++ )
    {
        
        if( data_tetrahedron_nextparents_of_edges[t][ei] != nullindex )
        assert( 0 <= data_tetrahedron_nextparents_of_edges[t][ei]
                &&
                data_tetrahedron_nextparents_of_edges[t][ei] < counter_tetrahedra );
    
        if( data_tetrahedron_nextparents_of_edges[t][ei] != nullindex )
        assert( data_tetrahedron_edges[ data_tetrahedron_nextparents_of_edges[t][ei] ][0] == data_tetrahedron_edges[t][ei] 
                ||
                data_tetrahedron_edges[ data_tetrahedron_nextparents_of_edges[t][ei] ][1] == data_tetrahedron_edges[t][ei] 
                ||
                data_tetrahedron_edges[ data_tetrahedron_nextparents_of_edges[t][ei] ][2] == data_tetrahedron_edges[t][ei] 
                ||
                data_tetrahedron_edges[ data_tetrahedron_nextparents_of_edges[t][ei] ][3] == data_tetrahedron_edges[t][ei] 
                ||
                data_tetrahedron_edges[ data_tetrahedron_nextparents_of_edges[t][ei] ][4] == data_tetrahedron_edges[t][ei] 
                ||
                data_tetrahedron_edges[ data_tetrahedron_nextparents_of_edges[t][ei] ][5] == data_tetrahedron_edges[t][ei] );
        
    }
    
    /* FACE TET STUFF */
    
    /* each tet: the next parents are unique           */
    /* each tet: the next parents are actually parents */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    for( int fi = 0; fi < 4; fi++ )
    {
        
        if( data_tetrahedron_nextparents_of_faces[t][fi] != nullindex )
            assert( 0 <= data_tetrahedron_nextparents_of_faces[t][fi] && data_tetrahedron_nextparents_of_faces[t][fi] < counter_faces );
    
        if( data_tetrahedron_nextparents_of_faces[t][fi] != nullindex )
            assert( data_tetrahedron_faces[ data_tetrahedron_nextparents_of_faces[t][fi] ][0] == data_tetrahedron_faces[t][fi] 
                    ||
                    data_tetrahedron_faces[ data_tetrahedron_nextparents_of_faces[t][fi] ][1] == data_tetrahedron_faces[t][fi] 
                    ||
                    data_tetrahedron_faces[ data_tetrahedron_nextparents_of_faces[t][fi] ][2] == data_tetrahedron_faces[t][fi] 
                    ||
                    data_tetrahedron_faces[ data_tetrahedron_nextparents_of_faces[t][fi] ][3] == data_tetrahedron_faces[t][fi] );
        
    }
    
    
    
    
    
    
   /****************************************/
   /* 6. Check that all parents are listed */
   /****************************************/
    
    
    
    /* check that each parent tetrahedron of a face is listed as a parent */
    
    for( int t  = 0; t  < counter_tetrahedra;  t++ )
    for( int fi = 0; fi <                  4; fi++ )
    {
      
      int f = data_tetrahedron_faces[t][fi];
      
      int p = data_face_firstparent_tetrahedron[f];
      
      assert( p != nullindex );
      
      while( p != t && p != nullindex )
        if( data_tetrahedron_faces[p][0] == f )
          p = data_tetrahedron_nextparents_of_faces[p][0];
        else if( data_tetrahedron_faces[p][1] == f )
          p = data_tetrahedron_nextparents_of_faces[p][1];
        else if( data_tetrahedron_faces[p][2] == f )
          p = data_tetrahedron_nextparents_of_faces[p][2];
        else if( data_tetrahedron_faces[p][3] == f )
          p = data_tetrahedron_nextparents_of_faces[p][3];
        else
          assert(false);
        
      assert( p == t );
      
    }
    
    /* check that each parent tetrahedron of an edge is listed as a parent */
    
    for( int t  = 0; t  < counter_tetrahedra; t++  )
    for( int ei = 0; ei <                  6; ei++ )
    {
      
      int e = data_tetrahedron_edges[t][ei];
      
      int p = data_edge_firstparent_tetrahedron[e];
      
      assert( p != nullindex );
      
      while( p != t && p != nullindex )
        if( data_tetrahedron_edges[p][0] == e )
          p = data_tetrahedron_nextparents_of_edges[p][0];
        else if( data_tetrahedron_edges[p][1] == e )
          p = data_tetrahedron_nextparents_of_edges[p][1];
        else if( data_tetrahedron_edges[p][2] == e )
          p = data_tetrahedron_nextparents_of_edges[p][2];
        else if( data_tetrahedron_edges[p][3] == e )
          p = data_tetrahedron_nextparents_of_edges[p][3];
        else if( data_tetrahedron_edges[p][4] == e )
          p = data_tetrahedron_nextparents_of_edges[p][4];
        else if( data_tetrahedron_edges[p][5] == e )
          p = data_tetrahedron_nextparents_of_edges[p][5];
        else
          assert(false);
        
      assert( p == t );
      
    }
    
    /* check that each parent tetrahedron of a vertex is listed as a parent */
    
    for( int t  = 0; t  < counter_tetrahedra;  t++ )
    for( int vi = 0; vi <                  4; vi++ )
    {
      
      int v = data_tetrahedron_vertices[t][vi];
      
      int p = data_vertex_firstparent_tetrahedron[v];
      
      assert( p != nullindex );
      
      while( p != t && p != nullindex )
        if( data_tetrahedron_vertices[p][0] == v )
          p = data_tetrahedron_nextparents_of_vertices[p][0];
        else if( data_tetrahedron_vertices[p][1] == v )
          p = data_tetrahedron_nextparents_of_vertices[p][1];
        else if( data_tetrahedron_vertices[p][2] == v )
          p = data_tetrahedron_nextparents_of_vertices[p][2];
        else if( data_tetrahedron_vertices[p][3] == v )
          p = data_tetrahedron_nextparents_of_vertices[p][3];
        else
          assert(false);
        
      assert( p == t );
      
    }
    
    /* check that each parent face of an edge is listed as a parent */
    
    for( int f  = 0; f  < counter_faces; f++ )
    for( int ei = 0; ei <             3; ei++ )
    {
      
      int e = data_face_edges[f][ei];
      
      int p = data_edge_firstparent_face[e];
      
      assert( p != nullindex );
      
      while( p != f && p != nullindex )
        if( data_face_edges[p][0] == e )
          p = data_face_nextparents_of_edges[p][0];
        else if( data_face_edges[p][1] == e )
          p = data_face_nextparents_of_edges[p][1];
        else if( data_face_edges[p][2] == e )
          p = data_face_nextparents_of_edges[p][2];
        else
          assert(false);
        
      assert( p == f );
      
    }
    
    /* check that each parent face of a vertex is listed as a parent */
    
    for( int f  = 0; f  < counter_faces;  f++ )
    for( int vi = 0; vi <             3; vi++ )
    {
      
      int v = data_face_vertices[f][vi];
      
      int p = data_vertex_firstparent_face[v];
      
      assert( p != nullindex );
      
      while( p != f && p != nullindex )
        if( data_face_vertices[p][0] == v )
          p = data_face_nextparents_of_vertices[p][0];
        else if( data_face_vertices[p][1] == v )
          p = data_face_nextparents_of_vertices[p][1];
        else if( data_face_vertices[p][2] == v )
          p = data_face_nextparents_of_vertices[p][2];
        else
          assert(false);
        
      assert( p == f );
      
    }
    
    /* check that each parent edge of a vertex is listed as a parent */
    
    for( int e  = 0; e  < counter_edges; e++ )
    for( int vi = 0; vi <             2; vi++ )
    {
      
      int v = data_edge_vertices[e][vi];
      
      int p = data_vertex_firstparent_edge[v];
      
      assert( p != nullindex );
      
      while( p != e && p != nullindex )
        p = data_edge_nextparents_of_vertices[p][ ( data_edge_vertices[p][0] == v ) ? 0 : 1 ];
      
      assert( p == e );
      
    }
    
    
    
    /*************************************/
    /* CHECK FROM BASE CLASS PERSPECTIVE */
    /*************************************/
   
    
    Mesh::check();
    
}
























void MeshSimplicial3D::print( std::ostream& os ) const
{
    os << "Printe Triangulation of 3D Manifold!" << std::endl;
    
    os << counter_tetrahedra << space << counter_faces << space << counter_edges << space << counter_vertices << nl;
    
    
    
    
    os << "Tetrahedron faces" << std::endl;
    
    for( const auto& quartett : data_tetrahedron_faces )
      std::cout << quartett[0] << space << quartett[1] << space << quartett[2] << space << quartett[3] << nl;
    
    os << "Faces first parent Tetrahedron" << std::endl;
    
    for( int fp : data_face_firstparent_tetrahedron )
      std::cout << fp << nl;
    
    os << "Tetrahedron next parents of faces" << std::endl;
    
    for( const auto& quartett : data_tetrahedron_nextparents_of_faces )
      std::cout << quartett[0] << space << quartett[1] << space << quartett[2] << space << quartett[3] << nl;
    
    
    
    os << "Tetrahedron edges" << std::endl;
    
    for( const auto& sextett : data_tetrahedron_edges )
      std::cout << sextett[0] << space << sextett[1] << space << sextett[2] << sextett[3] << space << sextett[4] << space << sextett[5] << nl;
    
    os << "Edge first parent tetrahedra" << std::endl;
    
    for( int fp : data_edge_firstparent_tetrahedron )
      std::cout << fp << nl;
    
    os << "Tetrahedron next parents of edges" << std::endl;
    
    for( const auto& sextett : data_tetrahedron_nextparents_of_edges )
      std::cout << sextett[0] << space << sextett[1] << space << sextett[2] << sextett[3] << space << sextett[4] << space << sextett[5] << nl;
    
    
    
    
    
    os << "Tetrahedron vertices" << std::endl;
    
    for( const auto& quartett : data_tetrahedron_vertices )
      std::cout << quartett[0] << space << quartett[1] << space << quartett[2] << space << quartett[3] << nl;
    
    os << "Edge first parent tetrahedra" << std::endl;
    
    for( int fp : data_vertex_firstparent_tetrahedron )
      std::cout << fp << nl;
    
    os << "Tetrahedron next parents of edges" << std::endl;
    
    for( const auto& quartett : data_tetrahedron_nextparents_of_vertices )
      std::cout << quartett[0] << space << quartett[1] << space << quartett[2] << space << quartett[3] << nl;
    
    
    
    
    
    
    os << "Face edges" << std::endl;
    
    for( const auto& triple : data_face_edges )
      std::cout << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    os << "Edge first parent faces" << std::endl;
    
    for( int fp : data_edge_firstparent_face )
      std::cout << fp << nl;
    
    os << "Face next parents of edges" << std::endl;
    
    for( const auto& triple : data_face_nextparents_of_edges )
      std::cout << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    
    
    
    
    os << "Face vertices" << std::endl;
    
    for( const auto& triple : data_face_vertices )
      std::cout << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    os << "Edge first parent faces" << std::endl;
    
    for( int fp : data_vertex_firstparent_face )
      std::cout << fp << nl;
    
    os << "Face next parents of edges" << std::endl;
    
    for( const auto& triple : data_face_nextparents_of_vertices )
      std::cout << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    
    
    
    os << "Edge vertices" << std::endl;
    
    for( const auto& duple : data_edge_vertices )
      std::cout << duple[0] << space << duple[1] << nl;
    
    os << "Vertex first parents" << std::endl;
    
    for( int fp : data_vertex_firstparent_edge )
      std::cout << fp << nl;
    
    os << "Edge next parents " << std::endl;
    
    for( const auto& duple : data_edge_nextparents_of_vertices )
      std::cout << duple[0] << space << duple[1] << nl;
    
    
    
    os << "Finished printing" << nl;
    
}









bool MeshSimplicial3D::dimension_counted( int dim ) const
{
    assert( 0 <= dim && dim <= 3 );
    return true;
}

int MeshSimplicial3D::count_simplices( int dim ) const
{
  if( dim == 0 )
    return count_vertices();
  else if( dim == 1 )
    return count_edges();
  else if( dim == 2 )
    return count_faces();
  else if( dim == 3 )
    return count_tetrahedra();
  else
    assert(false);
}

bool MeshSimplicial3D::subsimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub <= sup && sup <= 3 );
    return true;
}

IndexMap MeshSimplicial3D::getsubsimplices( int sup, int sub, int cell ) const
{
  
  if( sup == 3 && sub == 3 ) {
    
    assert( 0 <= cell && cell < count_tetrahedra() );
    return IndexMap( IndexRange(0,0), IndexRange(0,count_tetrahedra()-1), { cell } );
    
  } else if( sup == 3 && sub == 2 ) {
    
    assert( 0 <= cell && cell < count_tetrahedra() );
    auto temp = get_tetrahedron_faces(cell);
    return IndexMap( IndexRange(0,3), IndexRange(0,count_faces()-1), std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 3 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_tetrahedra() );
    auto temp = get_tetrahedron_edges(cell);
    return IndexMap( IndexRange(0,5), IndexRange( 0, count_edges()-1 ), std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 3 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_tetrahedra() );
    auto temp = get_tetrahedron_vertices(cell);
    return IndexMap( IndexRange(0,3), IndexRange( 0, count_vertices()-1 ) , std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 2 && sub == 2 ) {
    
    assert( 0 <= cell && cell < count_faces() );
    return IndexMap( IndexRange(0,0), IndexRange(0,count_faces()-1), { cell } );
    
  } else if( sup == 2 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_faces() );
    auto temp = get_face_edges(cell);
    return IndexMap( IndexRange(0,2), IndexRange( 0, count_edges()-1 ), std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 2 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_faces() );
    auto temp = get_face_vertices(cell);
    return IndexMap( IndexRange(0,2), IndexRange( 0, count_vertices()-1 ) , std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 1 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    return IndexMap( IndexRange(0,0), IndexRange(0,count_edges()-1), { cell } );
    
  } else if( sup == 1 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    auto temp = get_edge_vertices(cell);
    return IndexMap( IndexRange(0,1), IndexRange( 0, count_vertices()-1 ) , std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 0 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    return IndexMap( IndexRange(0,0), IndexRange(0,count_vertices()-1), { cell } );
    
  } else {
    
    assert(false);
    
  }
   
}

bool MeshSimplicial3D::supersimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub <= sup && sup <= 3 );
    return true;
}

const std::vector<int> MeshSimplicial3D::getsupersimplices( int sup, int sub, int cell ) const
{
  
  if( sup == 3 && sub == 3 ) {
    
    assert( 0 <= cell && cell < count_tetrahedra() );
    return { cell };
    
  } else if( sup == 3 && sub == 2 ) {
    
    assert( 0 <= cell && cell < count_faces() );
    auto temp = get_tetrahedron_parents_of_face( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 3 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    auto temp = get_tetrahedron_parents_of_edge( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 3 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    auto temp = get_tetrahedron_parents_of_vertex( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 2 && sub == 2 ) {
    
    assert( 0 <= cell && cell < count_faces() );
    return { cell };
    
  } else if( sup == 2 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    auto temp = get_face_parents_of_edge( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 2 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    auto temp = get_face_parents_of_vertex( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 1 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    return { cell };
    
  } else if( sup == 1 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    auto temp = get_edge_parents_of_vertex( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 0 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    return { cell };
    
  } else {
    
    assert(false);
    
  }
  
}








/* Count number of elements */

int MeshSimplicial3D::count_tetrahedra() const
{
    return counter_tetrahedra;
}

int MeshSimplicial3D::count_faces() const
{
    return counter_faces;
}

int MeshSimplicial3D::count_edges() const
{
    return counter_edges;
}

int MeshSimplicial3D::count_vertices() const
{
    return counter_vertices;
}



/* subsimplex relation of tetrahedra and faces */

bool MeshSimplicial3D::contains_tetrahedron_face( int t, int f ) const
{
    assert( 0 <= t && t < counter_tetrahedra );
    assert( 0 <= f && f < counter_faces );
    
    return ( data_tetrahedron_faces[t][0] == f ) 
           ||
           ( data_tetrahedron_faces[t][1] == f )
           ||
           ( data_tetrahedron_faces[t][2] == f ) 
           ||
           ( data_tetrahedron_faces[t][3] == f );
    
} 

int MeshSimplicial3D::indexof_tetrahedron_face( int t, int f ) const
{
    assert( 0 <= t && t < counter_tetrahedra );
    assert( 0 <= f && f < counter_faces );
    if     ( data_tetrahedron_faces[t][0] == f ) return 0;
    else if( data_tetrahedron_faces[t][1] == f ) return 1;
    else if( data_tetrahedron_faces[t][2] == f ) return 2;
    else if( data_tetrahedron_faces[t][3] == f ) return 3;
    else                                        assert(false);
} 

int MeshSimplicial3D::get_tetrahedron_face( int t, int fi ) const
{
    assert( 0 <= t  && t  < counter_tetrahedra );
    assert( 0 <= fi && fi < 4 );
    return data_tetrahedron_faces[t][fi];
} 

const std::array<int,4> MeshSimplicial3D::get_tetrahedron_faces( int t ) const
{
    assert( 0 <= t && t < counter_tetrahedra );
    return data_tetrahedron_faces[t];
} 



/* subsimplex relation of tetrahedra and edges */

bool MeshSimplicial3D::contains_tetrahedron_edge( int t, int e ) const
{
    assert( 0 <= t && t < counter_tetrahedra );
    assert( 0 <= e && e < counter_edges );
    
    return ( data_tetrahedron_edges[t][0] == e ) || ( data_tetrahedron_edges[t][1] == e ) || ( data_tetrahedron_edges[t][2] == e )
           ||
           ( data_tetrahedron_edges[t][3] == e ) || ( data_tetrahedron_edges[t][4] == e ) || ( data_tetrahedron_edges[t][5] == e )
           ;
} 

int MeshSimplicial3D::indexof_tetrahedron_edge( int t, int e ) const
{
    assert( 0 <= t && t < counter_tetrahedra );
    assert( 0 <= e && e < counter_edges );
    if     ( data_tetrahedron_edges[t][0] == e ) return 0;
    else if( data_tetrahedron_edges[t][1] == e ) return 1;
    else if( data_tetrahedron_edges[t][2] == e ) return 2;
    else if( data_tetrahedron_edges[t][3] == e ) return 3;
    else if( data_tetrahedron_edges[t][4] == e ) return 4;
    else if( data_tetrahedron_edges[t][5] == e ) return 5;
    else                                         assert(false);
} 

int MeshSimplicial3D::get_tetrahedron_edge( int t, int ei ) const
{
    assert( 0 <= t  && t  < counter_tetrahedra );
    assert( 0 <= ei && ei < 6 );
    return data_tetrahedron_edges[t][ei];
} 

const std::array<int,6> MeshSimplicial3D::get_tetrahedron_edges( int t ) const
{
    assert( 0 <= t && t < counter_tetrahedra );
    return data_tetrahedron_edges[t];
} 




/* subsimplex relation of tetrahedra and vertices */

bool MeshSimplicial3D::contains_tetrahedron_vertex( int t, int v ) const
{
    assert( 0 <= t && t < counter_tetrahedra );
    assert( 0 <= v && v < counter_vertices );
    
    return ( data_tetrahedron_vertices[t][0] == v ) 
           ||
           ( data_tetrahedron_vertices[t][1] == v )
           ||
           ( data_tetrahedron_vertices[t][2] == v )
           || 
           ( data_tetrahedron_vertices[t][3] == v );
    
} 

int MeshSimplicial3D::indexof_tetrahedron_vertex( int t, int v ) const
{
    assert( 0 <= t && t < counter_tetrahedra );
    assert( 0 <= v && v < counter_vertices );
    if     ( data_tetrahedron_vertices[t][0] == v ) return 0;
    else if( data_tetrahedron_vertices[t][1] == v ) return 1;
    else if( data_tetrahedron_vertices[t][2] == v ) return 2;
    else if( data_tetrahedron_vertices[t][3] == v ) return 3;
    else                                            assert(false);
} 

int MeshSimplicial3D::get_tetrahedron_vertex( int t, int vi ) const
{
    assert( 0 <= t  && t  < counter_tetrahedra );
    assert( 0 <= vi && vi < 4 );
    return data_tetrahedron_vertices[t][vi];
} 

const std::array<int,4> MeshSimplicial3D::get_tetrahedron_vertices( int t ) const
{
    assert( 0 <= t && t < counter_tetrahedra );
    return data_tetrahedron_vertices[t];
}



/* subsimplex relation of faces and edges */

bool MeshSimplicial3D::contains_face_edge( int f, int e ) const
{
    assert( 0 <= f && f < counter_faces );
    assert( 0 <= e && e < counter_edges );
    
    return ( data_face_edges[f][0] == e ) || ( data_face_edges[f][1] == e ) || ( data_face_edges[f][2] == e );
} 

int MeshSimplicial3D::indexof_face_edge( int f, int e ) const
{
    assert( 0 <= f && f < counter_faces );
    assert( 0 <= e && e < counter_edges );
    if     ( data_face_edges[f][0] == e ) return 0;
    else if( data_face_edges[f][1] == e ) return 1;
    else if( data_face_edges[f][2] == e ) return 2;
    else                                      assert(false);
} 

int MeshSimplicial3D::get_face_edge( int f, int ei ) const
{
    assert( 0 <= f  && f  < counter_faces );
    assert( 0 <= ei && ei < 3 );
    return data_face_edges[f][ei];
} 

const std::array<int,3> MeshSimplicial3D::get_face_edges( int f ) const
{
    assert( 0 <= f && f < counter_faces );
    return data_face_edges[f];
} 



/* subsimplex relation of face and vertices */

bool MeshSimplicial3D::contains_face_vertex( int f, int v ) const
{
    assert( 0 <= f && f < counter_faces );
    assert( 0 <= v && v < counter_vertices );
    
    return ( data_face_vertices[f][0] == v ) || ( data_face_vertices[f][1] == v ) || ( data_face_vertices[f][2] == v );
} 

int MeshSimplicial3D::indexof_face_vertex( int f, int v ) const
{
    assert( 0 <= f && f < counter_faces );
    assert( 0 <= v && v < counter_vertices );
    if     ( data_face_vertices[f][0] == v ) return 0;
    else if( data_face_vertices[f][1] == v ) return 1;
    else if( data_face_vertices[f][2] == v ) return 2;
    else                                     assert(false);
} 

int MeshSimplicial3D::get_face_vertex( int f, int vi ) const
{
    assert( 0 <= f  && f  < counter_faces );
    assert( 0 <= vi && vi < 3 );
    return data_face_vertices[f][vi];
} 

const std::array<int,3> MeshSimplicial3D::get_face_vertices( int f ) const
{
    assert( 0 <= f && f < counter_faces );
    return data_face_vertices[f];
} 




/* subsimplex relation of edges and vertices */

bool MeshSimplicial3D::contains_edge_vertex( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    
    return ( data_edge_vertices[e][0] == v ) || ( data_edge_vertices[e][1] == v );
} 

int MeshSimplicial3D::indexof_edge_vertex( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    if     ( data_edge_vertices[e][0] == v ) return 0;
    else if( data_edge_vertices[e][1] == v ) return 1;
    else                                     assert(false);
} 

int MeshSimplicial3D::get_edge_vertex( int e, int vi ) const
{
    assert( 0 <= e  && e  < counter_edges );
    assert( 0 <= vi && vi < 2 );
    return data_edge_vertices[e][vi];
} 

const std::array<int,2> MeshSimplicial3D::get_edge_vertices( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    return data_edge_vertices[e];
} 





/* tetrahedron parents of a face */

int MeshSimplicial3D::count_face_tetrahedron_parents( int f ) const
{
  return get_tetrahedron_parents_of_face( f ).size();
}

int MeshSimplicial3D::get_face_firstparent_tetrahedron( int f ) const
{
  assert( 0 <= f && f < counter_faces );
  return data_face_firstparent_tetrahedron[ f ];
}

int MeshSimplicial3D::get_face_nextparent_tetrahedron( int f, int t ) const
{
  assert( 0 <= f && f < counter_faces );
  assert( 0 <= t && t < counter_tetrahedra );
  
  if( data_tetrahedron_faces[t][0] == f )
    return data_tetrahedron_nextparents_of_faces[t][0];
  else if( data_tetrahedron_faces[t][1] == f )
    return data_tetrahedron_nextparents_of_faces[t][1];
  else if( data_tetrahedron_faces[t][2] == f )
    return data_tetrahedron_nextparents_of_faces[t][2];
  else if( data_tetrahedron_faces[t][3] == f )
    return data_tetrahedron_nextparents_of_faces[t][3];
  else
    assert(false);
}

int MeshSimplicial3D::get_tetrahedron_nextparent_of_face( int t, int fi ) const
{
  assert( 0 <= t  && t  < counter_tetrahedra );
  assert( 0 <= fi && fi < 4 );
  return data_tetrahedron_nextparents_of_faces[t][fi];
}

bool MeshSimplicial3D::is_tetrahedron_face_parent( int t, int f ) const
{
  assert( 0 <= f && f < counter_faces );
  assert( 0 <= t && t < counter_tetrahedra );
  return data_tetrahedron_faces[t][0] == f 
         ||
         data_tetrahedron_faces[t][1] == f
         ||
         data_tetrahedron_faces[t][2] == f
         ||
         data_tetrahedron_faces[t][3] == f;
}

int MeshSimplicial3D::indexof_tetrahedron_face_parent( int t, int f ) const
{
  assert( 0 <= f && f < count_faces() );
  std::vector<int> tetrahedra = get_tetrahedron_parents_of_face( f );
  
  auto iter = std::find( tetrahedra.begin(), tetrahedra.end(), t ); 
  assert( iter != tetrahedra.end() );
  
  return iter - tetrahedra.begin();
}

std::vector<int> MeshSimplicial3D::get_tetrahedron_parents_of_face( int f ) const
{
  assert( 0 <= f && f < count_faces() );
  
  std::vector<int> ret;
  
  for( int t = 0; t < count_tetrahedra(); t++ ) 
    for( int tf : get_tetrahedron_faces(t) )
      if( f == tf )
        ret.push_back( t );
  
  return ret;
}




/* tetrahedron parents of an edge */

int MeshSimplicial3D::count_edge_tetrahedron_parents( int e ) const
{
  return get_tetrahedron_parents_of_edge( e ).size();
}

int MeshSimplicial3D::get_edge_firstparent_tetrahedron( int e ) const
{
  assert( 0 <= e && e < counter_edges );
  return data_edge_firstparent_tetrahedron[ e ];
}

int MeshSimplicial3D::get_edge_nextparent_tetrahedron( int e, int t ) const
{
  assert( 0 <= e && e < counter_edges );
  assert( 0 <= t && t < counter_tetrahedra );
  
       if( data_tetrahedron_edges[t][0] == e )
    return data_tetrahedron_nextparents_of_edges[t][0];
  else if( data_tetrahedron_edges[t][1] == e )
    return data_tetrahedron_nextparents_of_edges[t][1];
  else if( data_tetrahedron_edges[t][2] == e )
    return data_tetrahedron_nextparents_of_edges[t][2];
  else if( data_tetrahedron_edges[t][3] == e )
    return data_tetrahedron_nextparents_of_edges[t][3];
  else if( data_tetrahedron_edges[t][4] == e )
    return data_tetrahedron_nextparents_of_edges[t][4];
  else if( data_tetrahedron_edges[t][5] == e )
    return data_tetrahedron_nextparents_of_edges[t][5];
  else
    assert(false);
}

int MeshSimplicial3D::get_tetrahedron_nextparent_of_edge( int t, int ei ) const
{
  assert( 0 <= t  && t  < counter_tetrahedra );
  assert( 0 <= ei && ei < 6 );
  return data_tetrahedron_nextparents_of_edges[t][ei];
}

bool MeshSimplicial3D::is_tetrahedron_edge_parent( int t, int e ) const
{
  assert( 0 <= e && e < counter_edges );
  assert( 0 <= t && t < counter_tetrahedra );
  return data_tetrahedron_edges[t][0] == e || data_tetrahedron_edges[t][1] == e || data_tetrahedron_edges[t][2] == e
         ||
         data_tetrahedron_edges[t][3] == e || data_tetrahedron_edges[t][4] == e || data_tetrahedron_edges[t][5] == e;
}

int MeshSimplicial3D::indexof_tetrahedron_edge_parent( int t, int e ) const
{
  assert( 0 <= e && e < count_edges() );
  std::vector<int> tetrahedra = get_tetrahedron_parents_of_edge( e );
  
  auto iter = std::find( tetrahedra.begin(), tetrahedra.end(), t ); 
  assert( iter != tetrahedra.end() );
  
  return iter - tetrahedra.begin();
}

std::vector<int> MeshSimplicial3D::get_tetrahedron_parents_of_edge( int e ) const
{
  assert( 0 <= e && e < count_edges() );
  
  std::vector<int> ret;
  
  for( int t = 0; t < count_tetrahedra(); t++ ) 
    for( int te : get_tetrahedron_edges(t) )
      if( e == te )
        ret.push_back( t );
  
  return ret;
}




/* tetrahedron parents of a vertex */

int MeshSimplicial3D::count_vertex_tetrahedron_parents( int v ) const
{
  return get_tetrahedron_parents_of_vertex( v ).size();
}

int MeshSimplicial3D::get_vertex_firstparent_tetrahedron( int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  return data_vertex_firstparent_tetrahedron[ v ];
}

int MeshSimplicial3D::get_vertex_nextparent_tetrahedron( int v, int t ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= t && t < counter_tetrahedra );
  
  if( data_tetrahedron_vertices[t][0] == v )
    return data_tetrahedron_nextparents_of_vertices[t][0];
  else if( data_tetrahedron_vertices[t][1] == v )
    return data_tetrahedron_nextparents_of_vertices[t][1];
  else if( data_tetrahedron_vertices[t][2] == v )
    return data_tetrahedron_nextparents_of_vertices[t][2];
  else if( data_tetrahedron_vertices[t][3] == v )
    return data_tetrahedron_nextparents_of_vertices[t][3];
  else
    assert(false);
}

int MeshSimplicial3D::get_tetrahedron_nextparent_of_vertex( int t, int vi ) const
{
  assert( 0 <= t  && t  < counter_tetrahedra );
  assert( 0 <= vi && vi < 4 );
  return data_tetrahedron_nextparents_of_vertices[t][vi];
}

bool MeshSimplicial3D::is_tetrahedron_vertex_parent( int t, int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= t && t < counter_tetrahedra );
  return data_tetrahedron_vertices[t][0] == v
         ||
         data_tetrahedron_vertices[t][1] == v
         ||
         data_tetrahedron_vertices[t][2] == v 
         ||
         data_tetrahedron_vertices[t][3] == v;
}

int MeshSimplicial3D::indexof_tetrahedron_vertex_parent( int t, int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  std::vector<int> tetrahedra = get_tetrahedron_parents_of_vertex( v );
  
  auto iter = std::find( tetrahedra.begin(), tetrahedra.end(), t ); 
  assert( iter != tetrahedra.end() );
  
  return iter - tetrahedra.begin();
}

std::vector<int> MeshSimplicial3D::get_tetrahedron_parents_of_vertex( int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  
  std::vector<int> ret;
  
  for( int t = 0; t < count_tetrahedra(); t++ ) 
    for( int tv : get_tetrahedron_vertices(t) )
      if( v == tv )
        ret.push_back( t );
  
  return ret;
}







/* face parents of a edge */

int MeshSimplicial3D::count_edge_face_parents( int e ) const
{
  return get_face_parents_of_edge( e ).size();
}

int MeshSimplicial3D::get_edge_firstparent_face( int e ) const
{
  assert( 0 <= e && e < counter_edges );
  return data_edge_firstparent_face[ e ];
}

int MeshSimplicial3D::get_edge_nextparent_face( int e, int f ) const
{
  assert( 0 <= e && e < counter_edges );
  assert( 0 <= f && f < counter_faces );
  
  if( data_face_edges[f][0] == e )
    return data_face_nextparents_of_edges[f][0];
  else if( data_face_edges[f][1] == e )
    return data_face_nextparents_of_edges[f][1];
  else if( data_face_edges[f][2] == e )
    return data_face_nextparents_of_edges[f][2];
  else
    assert(false);
}

int MeshSimplicial3D::get_face_nextparent_of_edge( int f, int ei ) const
{
  assert( 0 <= f  && f  < counter_faces );
  assert( 0 <= ei && ei < 3 );
  return data_face_nextparents_of_edges[f][ei];
}

bool MeshSimplicial3D::is_face_edge_parent( int f, int e ) const
{
  assert( 0 <= e && e < counter_edges );
  assert( 0 <= f && f < counter_faces );
  return data_face_edges[f][0] == e || data_face_edges[f][1] == e || data_face_edges[f][2] == e;
}

int MeshSimplicial3D::indexof_face_edge_parent( int f, int e ) const
{
  assert( 0 <= e && e < count_edges() );
  std::vector<int> faces = get_face_parents_of_edge( e );
  
  auto iter = std::find( faces.begin(), faces.end(), f ); 
  assert( iter != faces.end() );
  
  return iter - faces.begin();
}

std::vector<int> MeshSimplicial3D::get_face_parents_of_edge( int e ) const
{
  assert( 0 <= e && e < count_edges() );
  
  std::vector<int> ret;
  
  for( int f = 0; f < count_faces(); f++ ) 
    for( int fe : get_face_edges(f) )
      if( e == fe )
        ret.push_back( f );
  
  return ret;
}


/* face parents of a vertex */

int MeshSimplicial3D::count_vertex_face_parents( int v ) const
{
  return get_face_parents_of_vertex( v ).size();
}

int MeshSimplicial3D::get_vertex_firstparent_face( int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  return data_vertex_firstparent_face[ v ];
}

int MeshSimplicial3D::get_vertex_nextparent_face( int v, int f ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= f && f < counter_faces );
  
  if( data_face_vertices[f][0] == v )
    return data_face_nextparents_of_vertices[f][0];
  else if( data_face_vertices[f][1] == v )
    return data_face_nextparents_of_vertices[f][1];
  else if( data_face_vertices[f][2] == v )
    return data_face_nextparents_of_vertices[f][2];
  else
    assert(false);
}

int MeshSimplicial3D::get_face_nextparent_of_vertex( int f, int vi ) const
{
  assert( 0 <= f  && f  < counter_faces );
  assert( 0 <= vi && vi < 3 );
  return data_face_nextparents_of_vertices[f][vi];
}

bool MeshSimplicial3D::is_face_vertex_parent( int f, int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= f && f < counter_faces );
  return data_face_vertices[f][0] == v || data_face_vertices[f][1] == v || data_face_vertices[f][2] == v;
}

int MeshSimplicial3D::indexof_face_vertex_parent( int f, int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  std::vector<int> faces = get_face_parents_of_vertex( v );
  
  auto iter = std::find( faces.begin(), faces.end(), f ); 
  assert( iter != faces.end() );
  
  return iter - faces.begin();
}

std::vector<int> MeshSimplicial3D::get_face_parents_of_vertex( int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  
  std::vector<int> ret;
  
  for( int f = 0; f < count_faces(); f++ ) 
    for( int fv : get_face_vertices(f) )
      if( v == fv )
        ret.push_back( f );
  
  return ret;
}




/* edge parents of a vertex */

int MeshSimplicial3D::count_vertex_edge_parents( int v ) const
{
  return get_edge_parents_of_vertex( v ).size();
}

int MeshSimplicial3D::get_vertex_firstparent_edge( int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  return data_vertex_firstparent_edge[ v ];
}

int MeshSimplicial3D::get_vertex_nextparent_edge( int v, int e ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= e && e < counter_edges    );
  
  if( data_edge_vertices[e][0] == v )
    return data_edge_nextparents_of_vertices[e][0];
  else if( data_edge_vertices[e][1] == v )
    return data_edge_nextparents_of_vertices[e][1];
  else
    assert(false);
}

int MeshSimplicial3D::get_edge_nextparent_of_vertex( int e, int vi ) const
{
  assert( 0 <= e  && e  < counter_edges );
  assert( 0 <= vi && vi < 2 );
  return data_edge_nextparents_of_vertices[e][vi];
}

bool MeshSimplicial3D::is_edge_vertex_parent( int e, int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= e && e < counter_edges    );
  return data_edge_vertices[e][0] == v || data_edge_vertices[e][1] == v;
}

int MeshSimplicial3D::indexof_edge_vertex_parent( int e, int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  std::vector<int> edges = get_edge_parents_of_vertex( v );
  
  auto iter = std::find( edges.begin(), edges.end(), e ); 
  assert( iter != edges.end() );
  
  return iter - edges.begin();
}

std::vector<int> MeshSimplicial3D::get_edge_parents_of_vertex( int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  
  std::vector<int> ret;
  
  for( int e = 0; e < count_edges(); e++ ) 
    for( int ev : get_edge_vertices(e) )
      if( v == ev )
        ret.push_back( e );
  
  return ret;
}







FloatVector MeshSimplicial3D::get_tetrahedron_midpoint( int t ) const
{
    assert( 0 <= t && t < counter_tetrahedra );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = (   getcoordinates().getdata( get_tetrahedron_vertices(t)[0], d ) 
                 + getcoordinates().getdata( get_tetrahedron_vertices(t)[1], d ) 
                 + getcoordinates().getdata( get_tetrahedron_vertices(t)[2], d ) 
                 + getcoordinates().getdata( get_tetrahedron_vertices(t)[3], d )
                ) / 4.;
    return mid;
}

FloatVector MeshSimplicial3D::get_face_midpoint( int f ) const
{
    assert( 0 <= f && f < counter_faces );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = (   getcoordinates().getdata( get_face_vertices(f)[0], d ) 
                 + getcoordinates().getdata( get_face_vertices(f)[1], d ) 
                 + getcoordinates().getdata( get_face_vertices(f)[2], d )
                ) / 3.;
    return mid;
}

FloatVector MeshSimplicial3D::get_edge_midpoint    ( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = (   getcoordinates().getdata( get_edge_vertices(e)[0], d ) 
                 + getcoordinates().getdata( get_edge_vertices(e)[1], d )
                ) / 2.;
    return mid;
}







#include "mesh.simplicial3D.br.cpp"

#include "mesh.simplicial3D.ur.cpp"

#include "mesh.simplicial3D.mpr.cpp"






