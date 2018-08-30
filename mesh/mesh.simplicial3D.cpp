
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
    
    /* 1. Check the array sizes */ // OK
    
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
    
    
    
    /* VERTEX EDGE STUFF */
    
    /* 
     * each edge: each vertex is a valid index
     * each edge: each vertex is unique 
     * each edge: the next parents make sense 
     * each edge: the next parents are actually parents
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
        assert( data_edge_vertices[e][0] != nullindex );
        assert( data_edge_vertices[e][1] != nullindex );
        assert( 0 <= data_edge_vertices[e][0] && data_edge_vertices[e][0] < counter_vertices );
        assert( 0 <= data_edge_vertices[e][1] && data_edge_vertices[e][0] < counter_vertices );
        assert( data_edge_vertices[e][0] != data_edge_vertices[e][1] );
        
//         if( data_edge_nextparents_of_vertices[e][0] != nullindex || data_edge_nextparents_of_vertices[e][1] != nullindex )
//           assert( data_edge_nextparents_of_vertices[e][0] != data_edge_nextparents_of_vertices[e][1] );
        
          
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
        
    }
    
    
    
    
    
    
    /* VERTEX FACE STUFF */
    
    /*
     * each face: each vertex is a valid index
     * each face: each vertex is unique 
     * each face: the next parents are unique
     * each face: the next parents are actually parents 
     */
    
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
        
        
    }
    
    
    
    
    
    
    
    /* EDGE FACE STUFF */
    
    
    /* 
     * each face: each edge is a valid index
     * each face: each edge is unique 
     * each face: the next parents are unique
     * each face: the next parents are actually parents 
     */
    
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
        
        
        if( data_face_nextparents_of_edges[f][0] != nullindex || data_face_nextparents_of_edges[f][1] != nullindex )
          assert( data_face_nextparents_of_edges[f][0] != data_face_nextparents_of_edges[f][1] );
        
        if( data_face_nextparents_of_edges[f][0] != nullindex || data_face_nextparents_of_edges[f][2] != nullindex )
          assert( data_face_nextparents_of_edges[f][0] != data_face_nextparents_of_edges[f][2] );
        
        if( data_face_nextparents_of_edges[f][1] != nullindex || data_face_nextparents_of_edges[f][2] != nullindex )
          assert( data_face_nextparents_of_edges[f][1] != data_face_nextparents_of_edges[f][2] );
        
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
        
        
    }
    
    
    
    
    
    /* VERTEX TET STUFF */
    
    /* 
     * each tet: each vertex is a valid index
     * each tet: each vertex is unique 
     * each tet: the next parents are unique
     * each tet: the next parents are actually parents 
     */
    
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
        
        
    }
    
    
    
    /* EDGE TET STUFF */
    
    /* 
     * each tet: each edge is a valid index
     * each tet: each edge is unique 
     * each tet: the next parents are unique
     * each tet: the next parents are actually parents 
     */
    
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
        
        
    }
    
    
    
    
    
    /* FACE TET STUFF */
    
    /* 
     * each tet: each face is a valid index
     * each tet: each face is unique 
     * each tet: the next parents are unique
     * each tet: the next parents are actually parents 
     */
    
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
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
    
    
    
    
    
    /********************************************************/
    
    
    
    
    /* 
     * each first parent tetrahedron of an face: first parent is non-null and a valid parent
     */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        int p = data_face_firstparent_tetrahedron[t];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_tetrahedra );
        
        assert( data_tetrahedron_faces[p][0] == t || data_tetrahedron_faces[p][1] == t || data_tetrahedron_faces[p][2] == t || data_tetrahedron_faces[p][3] == t );
    }
    
    /* 
     * each first parent tetrahedron of an edge: first parent is non-null and a valid parent
     */
    
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
    
    /* 
     * each first parent tetrahedron of an vertex: first parent is non-null and a valid parent
     */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        int p = data_vertex_firstparent_tetrahedron[t];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_tetrahedra );
        
        assert( data_tetrahedron_vertices[p][0] == t || data_tetrahedron_vertices[p][1] == t || data_tetrahedron_vertices[p][2] == t || data_tetrahedron_vertices[p][3] == t );
    }
    
    /* 
     * each first parent face of an edge: first parent is non-null and a valid parent
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
        int p = data_edge_firstparent_face[e];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_faces );
        
        assert( data_face_edges[p][0] == e || data_face_edges[p][1] == e || data_face_edges[p][2] == e );
    }
    
    /* 
     * each first parent face of a vertex: first parent is non-null and a valid parent
     */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
        int p = data_vertex_firstparent_face[v];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_faces );
        
        assert( data_face_vertices[p][0] == v || data_face_vertices[p][1] == v || data_face_vertices[p][2] == v );
    }
    
    /* 
     * each first parent edge of a vertex: first parent is non-null and a valid parent
     */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
        int p = data_vertex_firstparent_edge[v];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_edges );
        
        assert( data_edge_vertices[p][0] == v || data_edge_vertices[p][1] == v );
    }
    
    
    /********************************************************/
    
    
    /* 
     * check that each parent tetrahedron of a face is listed as a parent 
     */
    
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
    
    /*
     * check that each parent tetrahedron of an edge is listed as a parent 
     */
    
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
    
    /* 
     * check that each parent tetrahedron of a vertex is listed as a parent 
     */
    
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
    
    /*
     * check that each parent face of an edge is listed as a parent 
     */
    
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
    
    /* 
     * check that each parent face of a vertex is listed as a parent 
     */
    
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
    
    /* 
     * check that each parent edge of a vertex is listed as a parent 
     */
    
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





/*
 * * * * * VERTEX BISECTION
 */
    
void MeshSimplicial3D::bisect_edge( int e )
{
    assert( 0 <= e && e < counter_edges );
    check();
    
    /*
     * DATA COLLECTION 
     */
    
    int e_back_vertex  = data_edge_vertices[ e ][ 0 ];
    int e_front_vertex = data_edge_vertices[ e ][ 1 ];
    
    std::vector<int> old_faces      = get_face_parents_of_edge( e );
    std::vector<int> old_tetrahedra = get_tetrahedron_parents_of_edge( e );
    
    std::vector<int> localindex_of_face_refinementedge( old_faces.size() );
    std::vector<int> localindex_of_tetrahedron_refinementedge( old_tetrahedra.size() );
    
    
    FloatVector midcoordinate = get_edge_midpoint( e );
    
    
    
    /*
     * ALLOCATE MEMORY FOR THE DATA  
     */
    
    data_tetrahedron_nextparents_of_faces.resize( counter_tetrahedra + old_tetrahedra.size(),               { nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_faces.resize               ( counter_tetrahedra + old_tetrahedra.size(),               { nullindex, nullindex, nullindex, nullindex } );
    data_face_firstparent_tetrahedron.resize    ( counter_faces + old_faces.size() + old_tetrahedra.size(),                                      nullindex );
    
    data_tetrahedron_nextparents_of_edges.resize( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_edges.resize               ( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    data_edge_firstparent_tetrahedron.resize    ( counter_edges + old_faces.size() + 1,                                                                  nullindex );
    
    data_tetrahedron_nextparents_of_vertices.resize( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_vertices.resize               ( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex } );
    data_vertex_firstparent_tetrahedron.resize     ( counter_vertices + 1,                                                            nullindex );
    
    data_face_nextparents_of_edges.resize( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_face_edges.resize               ( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_edge_firstparent_face.resize    ( counter_edges     + old_faces.size() + 1,                                           nullindex );
    
    data_face_nextparents_of_vertices.resize( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_face_vertices.resize               ( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_vertex_firstparent_face.resize     ( counter_vertices + 1,                                                               nullindex );
    
    data_edge_nextparents_of_vertices.resize( counter_edges    + old_faces.size() + 1, { nullindex, nullindex } );
    data_edge_vertices.resize               ( counter_edges    + old_faces.size() + 1, { nullindex, nullindex } );
    data_vertex_firstparent_edge.resize     ( counter_vertices + 1,                                   nullindex );
    
    
    
    
    /* TODO: Bis hierhin geupdated ^^^^^^ */
    /*
     * FILL IN DATA
     */
    
    /* vertices of the bisected edge */
    data_edge_vertices[ e ][ 0 ] = e_back_vertex;
    data_edge_vertices[ e ][ 1 ] = counter_vertices;
      
    data_edge_vertices[ counter_edges ][ 0 ] = counter_vertices;
    data_edge_vertices[ counter_edges ][ 1 ] = e_front_vertex;
    
    /* next parent of back vertex stays the same */
    /* next parent of front vertex */
    data_edge_nextparents_of_vertices[ counter_edges ][ 1 ] = data_edge_nextparents_of_vertices[ e ][ 1 ];
    
    // edge parent list of back vertex stays the same 
    // run over the front vertex edge parent list and replace 'e' by 'counter_edges'
    if( data_vertex_firstparent_edge[ e_front_vertex ] == e ) {
      
      data_vertex_firstparent_edge[ e_front_vertex ] = counter_edges;
      
    } else {
      
      int current_edge = data_vertex_firstparent_edge[ e_front_vertex ];
      
      while( data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] != e )
        current_edge = data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ];
      
      data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] = counter_edges;
      
    }
    
    /* first and next parents of new vertex */
    data_vertex_firstparent_edge[ counter_vertices ] = e;
    data_edge_nextparents_of_vertices[ e ][ 1 ] = counter_edges;
    data_edge_nextparents_of_vertices[ counter_edges ][ 0 ] = nullindex;
    
    // no parent faces yet for the new vertex; will be filled in below */
    data_vertex_firstparent_face[ counter_vertices ] = nullindex;
    
    // no parent faces yet for the front edge; will be filled in below */
    data_edge_firstparent_face[ counter_edges ] = nullindex;
    
    
    
    std::cout << "FIRST BIG LOOP" << std::endl;
    
    for( int of = 0; of < old_faces.size(); of++ ) {
      
      int f_old = old_faces[ of ];
      int f_new = counter_faces + of;
      
      int f_e0 = data_face_edges[ f_old ][ 0 ];
      int f_e1 = data_face_edges[ f_old ][ 1 ];
      int f_e2 = data_face_edges[ f_old ][ 2 ];
      
      int f_v0 = data_face_vertices[ f_old ][ 0 ];
      int f_v1 = data_face_vertices[ f_old ][ 1 ];
      int f_v2 = data_face_vertices[ f_old ][ 2 ];
      
      int f_e_n0 = data_face_nextparents_of_edges[ f_old ][ 0 ];
      int f_e_n1 = data_face_nextparents_of_edges[ f_old ][ 1 ];
      int f_e_n2 = data_face_nextparents_of_edges[ f_old ][ 2 ];
      
      int f_v_n0 = data_face_nextparents_of_vertices[ f_old ][ 0 ];
      int f_v_n1 = data_face_nextparents_of_vertices[ f_old ][ 1 ];
      int f_v_n2 = data_face_nextparents_of_vertices[ f_old ][ 2 ];
      
      localindex_of_face_refinementedge[ of ] = ( f_e0 == e ) ? 0 : ( ( f_e1 == e ) ? 1 : 2 );
      assert( data_face_edges[ f_old ][ localindex_of_face_refinementedge[ of ] ] == e );
      
      
      
      if( localindex_of_face_refinementedge[ of ] == 0 ) { // 0 1 
        
        assert( f_v0 == e_back_vertex && f_v1 == e_front_vertex && e == f_e0 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = f_v0;
        data_face_vertices[ f_old ][1] = counter_vertices;
        data_face_vertices[ f_old ][2] = f_v2;
        
        data_face_vertices[ f_new ][0] = counter_vertices;
        data_face_vertices[ f_new ][1] = f_v1;
        data_face_vertices[ f_new ][2] = f_v2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = e;
        data_face_edges[ f_old ][1] = f_e1;
        data_face_edges[ f_old ][2] = counter_edges + 1 + of;
        
        data_face_edges[ f_new ][0] = counter_edges;
        data_face_edges[ f_new ][1] = counter_edges + 1 + of;
        data_face_edges[ f_new ][2] = f_e2;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = counter_vertices;
        data_edge_vertices[ counter_edges + 1 + of ][1] = f_v2;
        
        
        /* face vertex neighbors */
        data_face_nextparents_of_vertices[ f_old ][0] = f_v_n0;
        data_face_nextparents_of_vertices[ f_old ][1] = counter_faces + of;
        data_face_nextparents_of_vertices[ f_old ][2] = counter_faces + of;
        
        data_face_nextparents_of_vertices[ f_new ][0] = data_vertex_firstparent_face[ counter_vertices ];
        data_face_nextparents_of_vertices[ f_new ][1] = f_v_n1;
        data_face_nextparents_of_vertices[ f_new ][2] = f_v_n2;
        
        data_vertex_firstparent_face[ counter_vertices ] = f_old;
        
        /* face edge neighbors */
        data_face_nextparents_of_edges[ f_old ][0] = f_e_n0;
        data_face_nextparents_of_edges[ f_old ][1] = f_e_n1;
        data_face_nextparents_of_edges[ f_old ][2] = counter_faces + of;
        
        data_face_nextparents_of_edges[ f_new ][0] = data_edge_firstparent_face[ counter_edges ];
        data_face_nextparents_of_edges[ f_new ][1] = nullindex;
        data_face_nextparents_of_edges[ f_new ][2] = f_e_n2;
        
        data_edge_firstparent_face[ counter_edges ] = f_new;
        
        data_edge_firstparent_face[ counter_edges + 1 + of ] = f_old;
        
        /* run over the front outer edge parent face list and replace 'f_old' by 'f_new' */
        if( data_edge_firstparent_face[ f_e2 ] == f_old ) 
          data_edge_firstparent_face[ f_e2 ] = f_new;
        else {
          int current_face = data_edge_firstparent_face[ f_e2 ];
          while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != f_old 
                 &&
                 data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != nullindex )
            current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ];
          assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != nullindex );
          data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] = f_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][0] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ f_v2 ];
        data_vertex_firstparent_edge[ f_v2 ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][1] = opposing_vertex_firstparent_edge;
        
        
        
      } else if( localindex_of_face_refinementedge[ of ] == 1 ) { // 0 2 
        
        assert( f_v0 == e_back_vertex && f_v2 == e_front_vertex && e == f_e1 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = f_v0;
        data_face_vertices[ f_old ][1] = f_v1;
        data_face_vertices[ f_old ][2] = counter_vertices;
        
        data_face_vertices[ f_new ][0] = f_v1;
        data_face_vertices[ f_new ][1] = counter_vertices;
        data_face_vertices[ f_new ][2] = f_v2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = f_e0;
        data_face_edges[ f_old ][1] = e;
        data_face_edges[ f_old ][2] = counter_edges + 1 + of;
        
        data_face_edges[ f_new ][0] = counter_edges + 1 + of;
        data_face_edges[ f_new ][1] = f_e2;
        data_face_edges[ f_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = f_v1;
        data_edge_vertices[ counter_edges + 1 + of ][1] = counter_vertices;
        
        
        /* face vertex neighbors */
        data_face_nextparents_of_vertices[ f_old ][0] = f_v_n0;
        data_face_nextparents_of_vertices[ f_old ][1] = counter_faces + of;
        data_face_nextparents_of_vertices[ f_old ][2] = counter_faces + of;
        
        data_face_nextparents_of_vertices[ f_new ][0] = f_v_n1;
        data_face_nextparents_of_vertices[ f_new ][1] = data_vertex_firstparent_face[ counter_vertices ];
        data_face_nextparents_of_vertices[ f_new ][2] = f_v_n2;
        
        data_vertex_firstparent_face[ counter_vertices ] = f_old;
        
        /* face edge neighbors */
        data_face_nextparents_of_edges[ f_old ][0] = f_e_n0;
        data_face_nextparents_of_edges[ f_old ][1] = f_e_n1;
        data_face_nextparents_of_edges[ f_old ][2] = counter_faces + of;
        
        data_face_nextparents_of_edges[ f_new ][0] = nullindex;
        data_face_nextparents_of_edges[ f_new ][1] = f_e_n2;
        data_face_nextparents_of_edges[ f_new ][2] = data_edge_firstparent_face[ counter_edges ];
        
        data_edge_firstparent_face[ counter_edges ] = f_new;
        
        data_edge_firstparent_face[ counter_edges + 1 + of ] = f_old;
        
        /* run over the front outer edge parent face list and replace 'f_old' by 'f_new' */
        if( data_edge_firstparent_face[ f_e2 ] == f_old ) 
          data_edge_firstparent_face[ f_e2 ] = f_new;
        else {
          int current_face = data_edge_firstparent_face[ f_e2 ];
          while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != f_old 
                 &&
                 data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != nullindex )
            current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ];
          assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != nullindex );
          data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] = f_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][1] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ f_v1 ];
        data_vertex_firstparent_edge[ f_v1 ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][0] = opposing_vertex_firstparent_edge;
        
        
      } else if( localindex_of_face_refinementedge[ of ] == 2 ) { // 1 2 
        
        assert( f_v1 == e_back_vertex && f_v2 == e_front_vertex && e == f_e2 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = f_v0;
        data_face_vertices[ f_old ][1] = f_v1;
        data_face_vertices[ f_old ][2] = counter_vertices;
        
        data_face_vertices[ f_new ][0] = f_v0;
        data_face_vertices[ f_new ][1] = counter_vertices;
        data_face_vertices[ f_new ][2] = f_v2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = f_e0;
        data_face_edges[ f_old ][1] = counter_edges + 1 + of;
        data_face_edges[ f_old ][2] = e;
        
        data_face_edges[ f_new ][0] = counter_edges + 1 + of;
        data_face_edges[ f_new ][1] = f_e1;
        data_face_edges[ f_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = f_v0;
        data_edge_vertices[ counter_edges + 1 + of ][1] = counter_vertices;
        
        
        /* face vertex neighbors */
        data_face_nextparents_of_vertices[ f_old ][0] = counter_faces + of;
        data_face_nextparents_of_vertices[ f_old ][1] = f_v_n1;
        data_face_nextparents_of_vertices[ f_old ][2] = counter_faces + of;
        
        data_face_nextparents_of_vertices[ f_new ][0] = f_v_n0;
        data_face_nextparents_of_vertices[ f_new ][1] = data_vertex_firstparent_face[ counter_vertices ];
        data_face_nextparents_of_vertices[ f_new ][2] = f_v_n2;
        
        data_vertex_firstparent_face[ counter_vertices ] = f_old;
        
        /* face edge neighbors */
        data_face_nextparents_of_edges[ f_old ][0] = f_e_n0;
        data_face_nextparents_of_edges[ f_old ][1] = counter_faces + of;
        data_face_nextparents_of_edges[ f_old ][2] = f_e_n2;
        
        data_face_nextparents_of_edges[ f_new ][0] = nullindex;
        data_face_nextparents_of_edges[ f_new ][1] = f_e_n1;
        data_face_nextparents_of_edges[ f_new ][2] = data_edge_firstparent_face[ counter_edges ];
        
        data_edge_firstparent_face[ counter_edges ] = f_new;
        
        data_edge_firstparent_face[ counter_edges + 1 + of ] = f_old;
        
        /* run over the front outer edge parent face list and replace 'f_old' by 'f_new' */
        if( data_edge_firstparent_face[ f_e1 ] == f_old ) 
          data_edge_firstparent_face[ f_e1 ] = f_new;
        else {
          int current_face = data_edge_firstparent_face[ f_e1 ];
          while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e1 ) ] != f_old 
                 &&
                 data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e1 ) ] != nullindex )
            current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e1 ) ];
          assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e1 ) ] != nullindex );
          data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e1 ) ] = f_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][1] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ f_v0 ];
        data_vertex_firstparent_edge[ f_v0 ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][0] = opposing_vertex_firstparent_edge;
        
      } else {
        
        assert(false);
        
      } 
      
    }
    
    
    std::cout << "FIRST BIG LOOP END" << std::endl;
    
    
    /* Run over the front vertex' parent faces and conduct manipulations */

    int* pointer_to_index = &data_vertex_firstparent_face[ e_front_vertex ];
    
    while( *pointer_to_index != nullindex ) { 
      
      std::vector<int>::iterator it = std::find( old_faces.begin(), old_faces.end(), *pointer_to_index );
      
      if( it != old_faces.end() ) {
        
        assert( *pointer_to_index == *it );
        assert( *it == old_faces[ it - old_faces.begin() ] );
        
        *pointer_to_index = counter_faces + ( it - old_faces.begin() );
        
        std::cout << "manipulate" << std::endl;
        
      } else std::cout << "keep" << std::endl;
      
      int localindex_of_front_vertex = nullindex;
      if( data_face_vertices[ *pointer_to_index ][ 0 ] == e_front_vertex ) localindex_of_front_vertex = 0;
      if( data_face_vertices[ *pointer_to_index ][ 1 ] == e_front_vertex ) localindex_of_front_vertex = 1;
      if( data_face_vertices[ *pointer_to_index ][ 2 ] == e_front_vertex ) localindex_of_front_vertex = 2;
      assert( localindex_of_front_vertex != nullindex );
      
      pointer_to_index = &( data_face_nextparents_of_vertices[ *pointer_to_index ][ localindex_of_front_vertex ] );
      
    }
    
    getcoordinates().append( midcoordinate );
    
    
    
    
    /* TODO: Fill in the data for the case distinctions below */
    
    std::cout << "SECOND BIG LOOP BEGIN" << std::endl;
    
    for( int ot = 0; ot < old_tetrahedra.size(); ot++ ) {
      
      int t_old = old_tetrahedra[ ot ];
      int t_new = counter_tetrahedra + ot;
      
      int f_new = counter_faces + old_faces.size() + ot;
      
      int t_f0 = data_tetrahedron_faces[ t_old ][ 0 ];
      int t_f1 = data_tetrahedron_faces[ t_old ][ 1 ];
      int t_f2 = data_tetrahedron_faces[ t_old ][ 2 ];
      int t_f3 = data_tetrahedron_faces[ t_old ][ 3 ];
      
      int t_e0 = data_tetrahedron_edges[ t_old ][ 0 ];
      int t_e1 = data_tetrahedron_edges[ t_old ][ 1 ];
      int t_e2 = data_tetrahedron_edges[ t_old ][ 2 ];
      int t_e3 = data_tetrahedron_edges[ t_old ][ 3 ];
      int t_e4 = data_tetrahedron_edges[ t_old ][ 4 ];
      int t_e5 = data_tetrahedron_edges[ t_old ][ 5 ];
      
      int t_v0 = data_tetrahedron_vertices[ t_old ][ 0 ];
      int t_v1 = data_tetrahedron_vertices[ t_old ][ 1 ];
      int t_v2 = data_tetrahedron_vertices[ t_old ][ 2 ];
      int t_v3 = data_tetrahedron_vertices[ t_old ][ 3 ];
      
      int t_f_n0 = data_tetrahedron_nextparents_of_faces[ t_old ][ 0 ];
      int t_f_n1 = data_tetrahedron_nextparents_of_faces[ t_old ][ 1 ];
      int t_f_n2 = data_tetrahedron_nextparents_of_faces[ t_old ][ 2 ];
      int t_f_n3 = data_tetrahedron_nextparents_of_faces[ t_old ][ 3 ];
      
      int t_e_n0 = data_tetrahedron_nextparents_of_edges[ t_old ][ 0 ];
      int t_e_n1 = data_tetrahedron_nextparents_of_edges[ t_old ][ 1 ];
      int t_e_n2 = data_tetrahedron_nextparents_of_edges[ t_old ][ 2 ];
      int t_e_n3 = data_tetrahedron_nextparents_of_edges[ t_old ][ 3 ];
      int t_e_n4 = data_tetrahedron_nextparents_of_edges[ t_old ][ 4 ];
      int t_e_n5 = data_tetrahedron_nextparents_of_edges[ t_old ][ 5 ];
      
      int t_v_n0 = data_tetrahedron_nextparents_of_vertices[ t_old ][ 0 ];
      int t_v_n1 = data_tetrahedron_nextparents_of_vertices[ t_old ][ 1 ];
      int t_v_n2 = data_tetrahedron_nextparents_of_vertices[ t_old ][ 2 ];
      int t_v_n3 = data_tetrahedron_nextparents_of_vertices[ t_old ][ 3 ];
      
      localindex_of_tetrahedron_refinementedge[ ot ] = nullindex;
      if( t_e0 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 0;
      if( t_e1 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 1;
      if( t_e2 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 2;
      if( t_e3 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 3;
      if( t_e4 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 4;
      if( t_e5 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 5;
      assert( localindex_of_tetrahedron_refinementedge[ ot ] != nullindex );
      assert( data_tetrahedron_edges[ t_old ][ localindex_of_tetrahedron_refinementedge[ ot ] ] == e );
      
      
      
      if(        localindex_of_face_refinementedge[ ot ] == 0 ) { // 0 1 
        
        assert( t_v0 == e_back_vertex && t_v1 == e_front_vertex && e == t_e0 );
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = ;
        data_tetrahedron_faces[ t_old ][1] = ;
        data_tetrahedron_faces[ t_old ][2] = ;
        data_tetrahedron_faces[ t_old ][3] = ;
        
        data_tetrahedron_faces[ t_new ][0] = ;
        data_tetrahedron_faces[ t_new ][1] = ;
        data_tetrahedron_faces[ t_new ][2] = ;
        data_tetrahedron_faces[ t_new ][3] = ;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = ;
        data_tetrahedron_edges[ t_old ][1] = ;
        data_tetrahedron_edges[ t_old ][2] = ;
        data_tetrahedron_edges[ t_old ][3] = ;
        data_tetrahedron_edges[ t_old ][4] = ;
        data_tetrahedron_edges[ t_old ][5] = ;
        
        data_tetrahedron_edges[ t_new ][0] = ;
        data_tetrahedron_edges[ t_new ][1] = ;
        data_tetrahedron_edges[ t_new ][2] = ;
        data_tetrahedron_edges[ t_new ][3] = ;
        data_tetrahedron_edges[ t_new ][4] = ;
        data_tetrahedron_edges[ t_new ][5] = ;
        
        /* Tetrahedron -> Vertices */
        data_triangles_vertices[ t_old ][0] = t_v0;
        data_triangles_vertices[ t_old ][1] = counter_vertices;
        data_triangles_vertices[ t_old ][2] = t_v2;
        data_triangles_vertices[ t_old ][3] = t_v3;
        
        data_triangles_vertices[ t_new ][0] = counter_vertices;
        data_triangles_vertices[ t_new ][1] = t_v1;
        data_triangles_vertices[ t_new ][2] = t_v2;
        data_triangles_vertices[ t_new ][3] = t_v3;
        
        /* Face -> Edges */
        data_face_edges[ t_old ][0] = ;
        data_face_edges[ t_old ][1] = ;
        data_face_edges[ t_old ][2] = ;
        
        data_face_edges[ t_new ][0] = ;
        data_face_edges[ t_new ][1] = ;
        data_face_edges[ t_new ][2] = ;
        
        /* Face -> Vertices */
        data_face_vertices[ t_old ][0] = ;
        data_face_vertices[ t_old ][1] = ;
        data_face_vertices[ t_old ][2] = ;
        
        data_face_vertices[ t_new ][0] = ;
        data_face_vertices[ t_new ][1] = ;
        data_face_vertices[ t_new ][2] = ;
        
        /* Neighbors: Tetrahedron -> Faces */
        /* Neighbors: Tetrahedron -> Edges */
        /* Neighbors: Tetrahedron -> Vertices */
        
        /* Neighbors: Face -> Edges */
        /* Neighbors: Face -> Vertices */
        
        /* run over the front outer face parent tetrahedra and replace 't_old' by 't_new' */
        if( data_face_firstparent_tetrahedron[ t_f? ] == t_old ) 
          data_face_firstparent_tetrahedron[ t_f? ] = t_new;
        else {
          int current_tetrahedron = data_face_firstparent_tetrahedron[ t_f? ];
          while( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f? ) ] != t_old 
                 &&
                 data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f? ) ] != nullindex )
            current_tetrahedron = data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_e? ) ];
          assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f? ) ] != nullindex );
          data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f? ) ] = t_new;
        }
        
        /* add new parent tetrahedron for the new vertex */
        int new_vertex_firstparent_tetrahedron = data_vertex_firstparent_tetrahedron[ counter_vertices ];
        data_vertex_firstparent_tetrahedron[ counter_vertices ] = ;
        data_edge_nextparents_of_vertices[  ][ ] = new_vertex_firstparent_tetrahedron;
        
        /* add new parent face for the opposing edge and its vertices */
        
        
        
      } else if( localindex_of_face_refinementedge[ ot ] == 1 ) { // 0 2 
        
        assert( t_v0 == e_back_vertex && t_v2 == e_front_vertex && e == t_e1 );
        
      } else if( localindex_of_face_refinementedge[ ot ] == 2 ) { // 0 3 
        
        assert( t_v0 == e_back_vertex && t_v3 == e_front_vertex && e == t_e2 );
    
      } else if( localindex_of_face_refinementedge[ ot ] == 3 ) { // 1 2 
        
        assert( t_v1 == e_back_vertex && t_v2 == e_front_vertex && e == t_e3 );
    
      } else if( localindex_of_face_refinementedge[ ot ] == 4 ) { // 1 3 
        
        assert( t_v1 == e_back_vertex && t_v3 == e_front_vertex && e == t_e4 );
    
      } else if( localindex_of_face_refinementedge[ ot ] == 5 ) { // 2 3 
        
        assert( t_v2 == e_back_vertex && t_v3 == e_front_vertex && e == t_e5 );
    
      } else {
        
        assert( false );
        
      }
       
      
    }
      
    std::cout << "SECOND BIG LOOP END" << std::endl;
    
    
    
    /* SOME FINISHINGS HERE, such as updating the adjacent tetrahedra for front vertex and front edges */
    
    
    
    
    
    
    std::cout << "FINISHED" << std::endl;
    
    /*
     *  UPDATE COUNTERS 
     */
    
    counter_faces += old_faces.size();
    counter_edges     += 1 + old_faces.size();
    counter_vertices  += 1;
    
    
    
    
    
    /* Done */
    
    check();
    
}

/*
 * * * * * VERTEX BISECTION
 */











/*
 *
 *  Split Tetrahedra: 8T 
 *  
 *    00 11 22 33
 *       00 01 02 03
 *       01 11 12 13
 *       02 12 22 23
 *       03 13 23 33
 *       
 *       01 02 03 13
 *       01 02 12 13 
 *       02 03 13 23
 *       02 12 13 23
 *
 *  Split Faces:    4F
 *    
 *    00 11 22
 *       00 01 02
 *       01 11 12
 *       02 12 22
 *       01 02 12
 *    00 11 33
 *       00 01 03
 *       01 11 13
 *       03 13 33
 *       01 03 13
 *    00 22 33
 *       00 02 03
 *       02 22 23
 *       03 23 33
 *       02 03 23
 *    11 22 33
 *       11 12 13
 *       12 22 23
 *       13 23 33
 *       12 13 23
 *       
 *  Completely New Faces: 8T
 *       
 *       01 02 03
 *       01 12 13
 *       02 12 23
 *       03 13 23
 *       
 *       01 02 13
 *       02 03 13
 *       02 12 13
 *       02 13 23
 *       
 * 
 *  Split Edges: 2 E 
 *    
 *    00 11
 *       00 01
 *       01 11
 *    00 22
 *       00 02
 *       02 22
 *    00 33
 *       00 03
 *       03 33
 *    11 22
 *       11 12
 *       12 22
 *    11 33
 *       11 13
 *       13 33
 *    22 33
 *       22 23
 *       23 33
 *    
 *  Completely New Edge: 3F + 1T
 *       
 *       
 *       01 02
 *       01 12
 *       02 12
 *       
 *       01 03
 *       01 13
 *       03 13
 *       
 *       02 03
 *       02 23
 *       03 23
 *       
 *       12 13
 *       12 23
 *       13 23
 *       
 * 
 * 
 *       02 13
 *      
 *
 *
 *
 */




// TODO: Debug the uniform refinement method 

void MeshSimplicial3D::uniformrefinement()
{
    check();
    
    
    int new_counter_tetrahedra = 8 * counter_tetrahedra;
    int new_counter_faces      = 4 * counter_faces + 8 * counter_tetrahedra;
    int new_counter_edges      = 2 * counter_edges + 3  * counter_faces + 1 * counter_tetrahedra;
    int new_counter_vertices   = 1 * counter_vertices + 1 * counter_edges;
    
    
    /* resize the arrays */
    
    /* tetrahedron -> face */
    
    data_tetrahedron_faces.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    data_face_firstparent_tetrahedron.resize( new_counter_faces, nullindex );
    
    data_tetrahedron_nextparents_of_faces.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    
    /* tetrahedron -> edge */
    
    data_tetrahedron_edges.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    
    data_edge_firstparent_tetrahedron.resize( new_counter_edges, nullindex );
    
    data_tetrahedron_nextparents_of_edges.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    
    
    /* tetrahedron -> vertex */
    
    data_tetrahedron_vertices.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    data_vertex_firstparent_tetrahedron.resize( new_counter_vertices, nullindex );
    
    data_tetrahedron_nextparents_of_vertices.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    
    /* face -> edge */
    
    data_face_edges.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    data_edge_firstparent_face.resize( new_counter_edges, nullindex );
    
    data_face_nextparents_of_edges.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    
    /* face -> vertex */
    
    data_face_vertices.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    data_vertex_firstparent_face.resize( new_counter_vertices, nullindex );
    
    data_face_nextparents_of_vertices.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    
    /* edge -> vertex */
    
    data_edge_vertices.resize( new_counter_edges, { nullindex, nullindex } );
    
    data_vertex_firstparent_edge.resize( new_counter_vertices, nullindex );
    
    data_edge_nextparents_of_vertices.resize( new_counter_edges, { nullindex, nullindex } );
    
    
    /* coordinates */
    
    getcoordinates().addcoordinates( counter_edges );
    
    
    
    
    
    
    /* 0. create the new coordinates and fill them up */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      getcoordinates().loadvector( counter_vertices + e, get_edge_midpoint( e ) );
    }
    
    
    
    
    
    /********************************/
    /***   VERTICES AND EDGES    ****/
    /********************************/
    
    
    /*** TREAT THE OLD VERTICES AND THEIR CONNECTION TO EDGES ****/
    
    /* 1. for each old vertex, set the new parent edge */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      int p = data_vertex_firstparent_edge[v];
      
      assert( p != nullindex && 0 <= p && p < counter_edges );
      
      int vi = data_edge_vertices[p][0] == v ? 0 : 1;
      
      assert( data_edge_vertices[p][0] == v || data_edge_vertices[p][1] == v );
      assert( data_edge_vertices[p][vi] == v );
      
      data_vertex_firstparent_edge[v] = vi * counter_edges + p;
    }
    
    
    /* 2. for each old edge, relocate the data of the old vertices' old parent edges */
    
    for( int e  = 0; e  < counter_edges;  e++ )
    for( int vi = 0; vi <=            1; vi++ )
    {
      int q = data_edge_nextparents_of_vertices[e][vi];
      
      int v = data_edge_vertices[e][vi];
      
      assert( v != nullindex && 0 <= v && v < counter_vertices );
      
      if( q == nullindex ) {
        
        data_edge_nextparents_of_vertices[ vi * counter_edges + e ][vi] = nullindex;
        
      } else {
        
        assert( 0 <= q && q < counter_edges );
        
        int vinp = data_edge_vertices[q][0] == v ? 0 : 1;
        
        assert( data_edge_vertices[q][0] == v || data_edge_vertices[q][1] == v );
        assert( data_edge_vertices[q][vinp] == v );
        
        data_edge_nextparents_of_vertices[ vi * counter_edges + e ][vi] = vinp * counter_edges + q;
      
      } 
      
    }
    
    
    /*** TREAT THE NEW VERTICES AND THEIR CONNECTION TO EDGES ****/
    
    /* 1.
     * for each new vertex (which is in the middle of an old edge), 
     * set the first and second parent edge from the old edge 
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      data_vertex_firstparent_edge[counter_vertices + e] = e;
      
      data_edge_nextparents_of_vertices[ 0 * counter_edges + e ][1] = counter_edges + e;
      data_edge_nextparents_of_vertices[ 1 * counter_edges + e ][0] = nullindex;
    }
    
    
    
    /* 2.
     * for each old edge, run over the adjacent faces 
     * and add the corresponding new edges to the list of 
     * parent edges of new vertex.
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int f = data_edge_firstparent_face[e];
      
      while( f != nullindex ) {
        
        int ei   = nullindex;
        int e_1  = nullindex; 
        int e_2  = nullindex;
        int vi_1 = nullindex;
        int vi_2 = nullindex;
        
        // [ 01 02 12 ] -> [ 01 02 ] [ 01 12 ] [ 02 12 ]
        
        if(        data_face_edges[f][0] == e ) {
          ei = 0; e_1 = 0; e_2 = 1; vi_1 = 0; vi_2 = 0; 
        } else if( data_face_edges[f][1] == e ) {
          ei = 1; e_1 = 0; e_2 = 2; vi_1 = 1; vi_2 = 0; 
        } else if( data_face_edges[f][2] == e ) {
          ei = 2; e_1 = 1; e_2 = 2; vi_1 = 1; vi_2 = 1; 
        } else
          assert(false);
        
        assert( ei  != nullindex && e_1 != nullindex && e_2 != nullindex );
        
        int old_first_parent = data_vertex_firstparent_edge[ counter_vertices + e ];
        
        data_vertex_firstparent_edge[ counter_vertices + e ]
          = 2 * counter_edges + e_1 * counter_faces + f;
        
        data_edge_nextparents_of_vertices[ 2 * counter_edges + e_1 * counter_faces + f ][ vi_1 ]
          = 2 * counter_edges + e_2 * counter_faces + f;
        
        data_edge_nextparents_of_vertices[ 2 * counter_edges + e_2 * counter_faces + f ][ vi_2 ]
          = old_first_parent;
        
        f = data_face_nextparents_of_edges[ f ][ ei ];
        
      }
      
    }
    
    
    
    
    /* 3.
     * for each tetrahedron, include the single interior edge 02 13
     * in the parent lists of its two vertices
     */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      /* 01 02 03 12 13 23 */
      
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      
      int ne = 2 * counter_edges + 3 * counter_faces + t;
      
      int fp_v02 = data_vertex_firstparent_edge[ v02 ];
      int fp_v13 = data_vertex_firstparent_edge[ v13 ];
      
      data_vertex_firstparent_edge[ v02 ] = ne;
      data_vertex_firstparent_edge[ v13 ] = ne;
      
      data_edge_nextparents_of_vertices[ ne ][ 0 ] = fp_v02;
      data_edge_nextparents_of_vertices[ ne ][ 1 ] = fp_v13;
      
    }
    
    
    
    
    /*** SET THE VERTICES OF ALL EDGES ****/
    
    /* 1. for each edge created from an old edge, set the vertices */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int vertex_back  = data_edge_vertices[e][0];
      int vertex_front = data_edge_vertices[e][1];
      
      data_edge_vertices[e                ][0] = vertex_back;
      data_edge_vertices[e                ][1] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][0] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][1] = vertex_front;
    }
    
    /* 2. for each face, set the vertices of the new edges inside the face */
    
    for( int f = 0; f < counter_faces; f++ )
    {
      data_edge_vertices[ 2 * counter_edges + 0 * counter_faces + f ][0] = counter_vertices + data_face_edges[f][0];
      data_edge_vertices[ 2 * counter_edges + 0 * counter_faces + f ][1] = counter_vertices + data_face_edges[f][1];
      data_edge_vertices[ 2 * counter_edges + 1 * counter_faces + f ][0] = counter_vertices + data_face_edges[f][0];
      data_edge_vertices[ 2 * counter_edges + 1 * counter_faces + f ][1] = counter_vertices + data_face_edges[f][2];
      data_edge_vertices[ 2 * counter_edges + 2 * counter_faces + f ][0] = counter_vertices + data_face_edges[f][1];
      data_edge_vertices[ 2 * counter_edges + 2 * counter_faces + f ][1] = counter_vertices + data_face_edges[f][2];
    }
    
    /* 3. for each tetrahedron, set the internal vertices */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      /* 01 02 03 12 13 23 */
      
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      int ne = 2 * counter_edges + 3 * counter_faces + t;
      
      data_edge_vertices[ ne ][ 0 ] = v02;
      data_edge_vertices[ ne ][ 1 ] = v13;
      
    }
    
    
    
    
    
    /****************************************/
    /***   VERTICES AND EDGES AND FACES  ****/
    /****************************************/
    
    
    
    /*** TREAT THE OLD VERTICES AND THEIR CONNECTION TO FACES ****/
    
    /* 1. for each old vertex, set the new parent face */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      int p = data_vertex_firstparent_face[v];
      
      assert( p != nullindex && 0 <= p && p < counter_faces );
      
      int vi = data_face_vertices[p][0] == v ? 0 : data_face_vertices[p][1] == v ? 1 : 2;
      
      assert( data_face_vertices[p][0] == v || data_face_vertices[p][1] == v || data_face_vertices[p][2] == v );
      assert( data_face_vertices[p][vi] == v );
      
      data_vertex_firstparent_face[v] = vi * counter_faces + p;
    }
    
    
    /* 2. for each old face, relocate the data of the old vertices' parent face */
    
    for( int f  = 0; f  < counter_faces;  f++ )
    for( int vi = 0; vi <             3; vi++ )
    {
      
      int q = data_face_nextparents_of_vertices[f][vi];
      
      int v = data_face_vertices[f][vi];
      
      if( q == nullindex ) {
        
        data_face_nextparents_of_vertices[ vi * counter_faces + f ][vi] = nullindex;
        
      } else {
        
        int vinp = data_face_vertices[q][0] == v ? 0 : data_face_vertices[q][1] == v ? 1 : 2;
        
        assert( data_face_vertices[q][0] == v || data_face_vertices[q][1] == v || data_face_vertices[q][2] == v );
        assert( data_face_vertices[q][vinp] == v );
        
        data_face_nextparents_of_vertices[ vi * counter_faces + f ][vi] = vinp * counter_faces + q;
      
      } 
      
    }
    
    
    /*** TREAT THE NEW VERTICES AND THEIR CONNECTION TO FACES ****/
    
    /* 1.
     * for each old edge, run over the adjacent old faces 
     * and add the corresponding new faces to the list of 
     * parent faces of new vertex.
     */
   
   /* TODO: the following code is currently considered legacy */
    
//     for( int e = 0; e < counter_edges; e++ )
//     {
//       int f = data_edge_firstparent_face[e];
//       
//       while( f != nullindex ) {
//         
//         int ei   = nullindex;
//         int f_1  = nullindex; 
//         int f_2  = nullindex;
//         int f_3  = nullindex;
//         int vi_1 = nullindex;
//         int vi_2 = nullindex;
//         int vi_3 = nullindex;
//         
//         // [ 01 02 12 ] -> [ 01 02 ] [ 01 12 ] [ 02 12 ]
//         // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
//         
//         if(        data_face_edges[f][0] == e ) {
//           ei = 0; 
//           f_1 = 0; f_2 = 3; f_3 = 1; vi_1 = 1; vi_2 = 0; vi_3 = 0;  
//         } else if( data_face_edges[f][1] == e ) {
//           ei = 1; 
//           f_1 = 0; f_2 = 3; f_3 = 2; vi_1 = 2; vi_2 = 1; vi_3 = 0;  
//         } else if( data_face_edges[f][2] == e ) {
//           ei = 2; 
//           f_1 = 1; f_2 = 3; f_3 = 2; vi_1 = 2; vi_2 = 2; vi_3 = 1; 
//         } else
//           assert(false);
//         
//         int old_first_parent = data_vertex_firstparent_face[ counter_vertices + e ];
//         
//         data_vertex_firstparent_face[ counter_vertices + e ]
//           = f_1 * counter_faces + f;
//         
//         if( f_1 != 0 ) assert( data_face_nextparents_of_vertices[ f_1 * counter_faces + f ][ vi_1 ] == nullindex );
//         data_face_nextparents_of_vertices[ f_1 * counter_faces + f ][ vi_1 ]
//           = f_2 * counter_faces + f;
//         
//         if( f_2 != 0 ) assert( data_face_nextparents_of_vertices[ f_2 * counter_faces + f ][ vi_2 ] == nullindex );
//         data_face_nextparents_of_vertices[ f_2 * counter_faces + f ][ vi_2 ]
//           = f_3 * counter_faces + f;
//         
//         if( f_3 != 0 ) assert( data_face_nextparents_of_vertices[ f_3 * counter_faces + f ][ vi_3 ] == nullindex );
//         data_face_nextparents_of_vertices[ f_3 * counter_faces + f ][ vi_3 ]
//           = old_first_parent;
//         
//         f = data_face_nextparents_of_edges[ f ][ ei ];
//         
//       }
//       
//     }
    
    for( int f = 0; f < counter_faces; f++ )
    {
        
        int e01 = data_face_edges[ f ][0];
        int e02 = data_face_edges[ f ][1];
        int e12 = data_face_edges[ f ][2];
        
        int v01 = counter_vertices + e01;
        int v02 = counter_vertices + e02;
        int v12 = counter_vertices + e12;
        
        int f_00_01_02 = 0 * counter_faces + f;
        int f_01_11_12 = 1 * counter_faces + f;
        int f_02_12_22 = 2 * counter_faces + f;
        int f_01_02_12 = 3 * counter_faces + f;
        
        int fp_v01 = data_vertex_firstparent_face[ v01 ];
        int fp_v02 = data_vertex_firstparent_face[ v02 ];
        int fp_v12 = data_vertex_firstparent_face[ v12 ];
        
        // [ 01 02 12 ] -> [ 01 02 ] [ 01 12 ] [ 02 12 ]
        // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
        data_vertex_firstparent_face[ v01 ] = f_00_01_02;
        data_face_nextparents_of_vertices[ f_00_01_02 ][ 1 ] = f_01_11_12;
        data_face_nextparents_of_vertices[ f_01_11_12 ][ 0 ] = f_01_02_12;
        data_face_nextparents_of_vertices[ f_01_02_12 ][ 0 ] = fp_v01;
        
        data_vertex_firstparent_face[ v02 ] = f_00_01_02;
        data_face_nextparents_of_vertices[ f_00_01_02 ][ 2 ] = f_02_12_22;
        data_face_nextparents_of_vertices[ f_02_12_22 ][ 0 ] = f_01_02_12;
        data_face_nextparents_of_vertices[ f_01_02_12 ][ 1 ] = fp_v02;
        
        data_vertex_firstparent_face[ v12 ] = f_01_11_12;
        data_face_nextparents_of_vertices[ f_01_11_12 ][ 2 ] = f_02_12_22;
        data_face_nextparents_of_vertices[ f_02_12_22 ][ 1 ] = f_01_02_12;
        data_face_nextparents_of_vertices[ f_01_02_12 ][ 2 ] = fp_v12;
        
    }
    
    
    
    
    
    
    /* for each tetrahedron, add the new interior faces to the parent lists of the new vertices */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
       /*
        *       01 02 03
        *       01 12 13
        *       02 12 23
        *       03 13 23
        *       
        *       01 02 13
        *       02 03 13
        *       02 12 13
        *       02 13 23
        */
      
      int v01 = counter_vertices + data_tetrahedron_edges[t][0];
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v03 = counter_vertices + data_tetrahedron_edges[t][2];
      int v12 = counter_vertices + data_tetrahedron_edges[t][3];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      int v23 = counter_vertices + data_tetrahedron_edges[t][5];
      
      int fp_v01 = data_vertex_firstparent_face[ v01 ];
      int fp_v02 = data_vertex_firstparent_face[ v02 ];
      int fp_v03 = data_vertex_firstparent_face[ v03 ];
      int fp_v12 = data_vertex_firstparent_face[ v12 ];
      int fp_v13 = data_vertex_firstparent_face[ v13 ];
      int fp_v23 = data_vertex_firstparent_face[ v23 ];
      
      
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      
      data_vertex_firstparent_face[ v01 ] = f_01_02_03;
      data_face_nextparents_of_vertices[ f_01_02_03 ][ 0 ] = f_01_12_13; 
      data_face_nextparents_of_vertices[ f_01_12_13 ][ 0 ] = f_01_02_13; 
      data_face_nextparents_of_vertices[ f_01_02_13 ][ 0 ] = fp_v01; 
      
      data_vertex_firstparent_face[ v02 ] = f_01_02_03;
      data_face_nextparents_of_vertices[ f_01_02_03 ][ 1 ] = f_02_12_23; 
      data_face_nextparents_of_vertices[ f_02_12_23 ][ 0 ] = f_01_02_13; 
      data_face_nextparents_of_vertices[ f_01_02_13 ][ 1 ] = f_02_03_13; 
      data_face_nextparents_of_vertices[ f_02_03_13 ][ 0 ] = f_02_12_13; 
      data_face_nextparents_of_vertices[ f_02_12_13 ][ 0 ] = f_02_13_23; 
      data_face_nextparents_of_vertices[ f_02_13_23 ][ 0 ] = fp_v02; 
      
      data_vertex_firstparent_face[ v03 ] = f_01_02_03;
      data_face_nextparents_of_vertices[ f_01_02_03 ][ 2 ] = f_03_13_23; 
      data_face_nextparents_of_vertices[ f_03_13_23 ][ 0 ] = f_02_03_13; 
      data_face_nextparents_of_vertices[ f_02_03_13 ][ 1 ] = fp_v03; 
      
      data_vertex_firstparent_face[ v12 ] = f_01_12_13;
      data_face_nextparents_of_vertices[ f_01_12_13 ][ 1 ] = f_02_12_23; 
      data_face_nextparents_of_vertices[ f_02_12_23 ][ 1 ] = f_02_12_13; 
      data_face_nextparents_of_vertices[ f_02_12_13 ][ 1 ] = fp_v12; 
      
      data_vertex_firstparent_face[ v13 ] = f_01_12_13;
      data_face_nextparents_of_vertices[ f_01_12_13 ][ 2 ] = f_03_13_23; 
      data_face_nextparents_of_vertices[ f_03_13_23 ][ 1 ] = f_02_03_13; 
      data_face_nextparents_of_vertices[ f_02_03_13 ][ 2 ] = f_02_12_13; 
      data_face_nextparents_of_vertices[ f_02_12_13 ][ 2 ] = f_02_13_23; 
      data_face_nextparents_of_vertices[ f_02_13_23 ][ 1 ] = fp_v13; 
      
      data_vertex_firstparent_face[ v23 ] = f_02_12_23;
      data_face_nextparents_of_vertices[ f_02_12_23 ][ 2 ] = f_03_13_23; 
      data_face_nextparents_of_vertices[ f_03_13_23 ][ 2 ] = f_02_13_23; 
      data_face_nextparents_of_vertices[ f_02_13_23 ][ 2 ] = fp_v23; 
      
    }
    
    
    
    /*** TREAT THE BISECTED EDGES AND THEIR CONNECTION TO FACES ****/
    
    
    /* for each bisected edge, set the new first parent faces of the children edges */
    for( int e = 0; e < counter_edges; e++ )
    {
      int p = data_edge_firstparent_face[e];
      
      assert( p != nullindex );
      
      int ei        = nullindex;
      int nfp_back  = nullindex;
      int nfp_front = nullindex;
      
      if( data_face_edges[p][0] == e ){
        ei = 0; nfp_back = 0; nfp_front = 1;
      } else if( data_face_edges[p][1] == e ) {
        ei = 1; nfp_back = 0; nfp_front = 2;
      } else if( data_face_edges[p][2] == e ) {
        ei = 2; nfp_back = 1; nfp_front = 2;
      } else 
        assert(false);
      
      assert( ei != nullindex );
      assert( data_face_edges[p][0] == e || data_face_edges[p][1] == e || data_face_edges[p][2] == e );
      assert( data_face_edges[p][ei] == e );
      
      data_edge_firstparent_face[ 0 * counter_edges + e ] = nfp_back  * counter_faces + p;
      data_edge_firstparent_face[ 1 * counter_edges + e ] = nfp_front * counter_faces + p;
      
    }
    
    
    /* for each face, relocate the data of the old edges' parent face */
    /* additionally, set the new parents */
    for( int f  = 0; f  < counter_faces;  f++ )
    for( int ei = 0; ei <             3; ei++ )
    {
      int e = data_face_edges[f][ei];
      
      int f_back  = nullindex;
      int f_front = nullindex;
      int e_back  = nullindex;
      int e_front = nullindex;
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
      
      if( ei == 0 ){
        f_back = 0; f_front = 1; e_back = 0; e_front = 0; 
      } else if( ei == 1 ) {
        f_back = 0; f_front = 2; e_back = 1; e_front = 1; 
      } else if( ei == 2 ) {
        f_back = 1; f_front = 2; e_back = 2; e_front = 2; 
      } else 
        assert(false);
      
      
      int q = data_face_nextparents_of_edges[f][ei];
      
      if( q == nullindex ) {
        
        data_face_nextparents_of_edges[ f_back  * counter_faces + f ][ e_back  ] = nullindex;
        data_face_nextparents_of_edges[ f_front * counter_faces + f ][ e_front ] = nullindex;
        
      } else if( q != nullindex ) {
        
        int q_ei        = nullindex;
        int q_nfp_back  = nullindex;
        int q_nfp_front = nullindex;
        
        if(        data_face_edges[q][0] == e ){
          q_ei = 0; q_nfp_back = 0; q_nfp_front = 1;
        } else if( data_face_edges[q][1] == e ) {
          q_ei = 1; q_nfp_back = 0; q_nfp_front = 2;
        } else if( data_face_edges[q][2] == e ) {
          q_ei = 2; q_nfp_back = 1; q_nfp_front = 2;
        } else 
          assert(false);
        
        assert( q_ei != nullindex );
        assert( data_face_edges[q][0] == e || data_face_edges[q][1] == e || data_face_edges[q][2] == e );
        assert( data_face_edges[q][q_ei] == e );
        
        data_face_nextparents_of_edges[ f_back  * counter_faces + f ][ e_back  ] = q_nfp_back  * counter_faces + q;
        data_face_nextparents_of_edges[ f_front * counter_faces + f ][ e_front ] = q_nfp_front * counter_faces + q;
        
      } 
      
    }
    
    
    /*** TREAT THE FACE-BASED EDGES AND THEIR CONNECTION TO SPLIT FACES ****/
    
    /* for each face, run over the new edges and add firstparents and parents */
    for( int f = 0; f < counter_faces; f++ )
    {
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
      // new edges [ 01 02 ], [ 01 12 ], [ 02 12 ]
      
      data_edge_firstparent_face[ 2 * counter_edges + 0 * counter_faces + f ] = 3 * counter_faces + f;
      data_edge_firstparent_face[ 2 * counter_edges + 1 * counter_faces + f ] = 3 * counter_faces + f;
      data_edge_firstparent_face[ 2 * counter_edges + 2 * counter_faces + f ] = 3 * counter_faces + f;
      
      data_face_nextparents_of_edges[ 3 * counter_faces + f ][0] = 0 * counter_faces + f;
      data_face_nextparents_of_edges[ 3 * counter_faces + f ][1] = 1 * counter_faces + f;
      data_face_nextparents_of_edges[ 3 * counter_faces + f ][2] = 2 * counter_faces + f;
      
      data_face_nextparents_of_edges[ 0 * counter_faces + f ][2] = nullindex;
      data_face_nextparents_of_edges[ 1 * counter_faces + f ][1] = nullindex;
      data_face_nextparents_of_edges[ 2 * counter_faces + f ][0] = nullindex;
      
    }
    
    /*** TREAT THE FACE-BASED EDGES AND THEIR CONNECTION TO INTERIOR FACES ****/
    /*** TREAT THE SINGLE INTERIOR EDGE AND ITS CONNECTION TO INTERIOR FACES ****/
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
       /*
        *       01 02 03
        *       01 12 13
        *       02 12 23
        *       03 13 23
        *       
        *       01 02 13
        *       02 03 13
        *       02 12 13
        *       02 13 23
        */
      
      int f_012 = data_tetrahedron_faces[t][0];
      int f_013 = data_tetrahedron_faces[t][1];
      int f_023 = data_tetrahedron_faces[t][2];
      int f_123 = data_tetrahedron_faces[t][3];
      
      int e_01_02 = 2 * counter_edges + 0 * counter_faces + f_012;
      int e_01_12 = 2 * counter_edges + 1 * counter_faces + f_012;
      int e_02_12 = 2 * counter_edges + 2 * counter_faces + f_012;
      
      int e_01_03 = 2 * counter_edges + 0 * counter_faces + f_013;
      int e_01_13 = 2 * counter_edges + 1 * counter_faces + f_013;
      int e_03_13 = 2 * counter_edges + 2 * counter_faces + f_013;
      
      int e_02_03 = 2 * counter_edges + 0 * counter_faces + f_023;
      int e_02_23 = 2 * counter_edges + 1 * counter_faces + f_023;
      int e_03_23 = 2 * counter_edges + 2 * counter_faces + f_023;
      
      int e_12_13 = 2 * counter_edges + 0 * counter_faces + f_123;
      int e_12_23 = 2 * counter_edges + 1 * counter_faces + f_123;
      int e_13_23 = 2 * counter_edges + 2 * counter_faces + f_123;
      
      int e_02_13 = 2 * counter_edges + 3 * counter_faces + t;
      
      
      int fp_e_01_02 = data_edge_firstparent_face[ e_01_02 ];
      int fp_e_01_12 = data_edge_firstparent_face[ e_01_12 ];
      int fp_e_02_12 = data_edge_firstparent_face[ e_02_12 ];
      
      int fp_e_01_03 = data_edge_firstparent_face[ e_01_03 ];
      int fp_e_01_13 = data_edge_firstparent_face[ e_01_13 ];
      int fp_e_03_13 = data_edge_firstparent_face[ e_03_13 ];
      
      int fp_e_02_03 = data_edge_firstparent_face[ e_02_03 ];
      int fp_e_02_23 = data_edge_firstparent_face[ e_02_23 ];
      int fp_e_03_23 = data_edge_firstparent_face[ e_03_23 ];
      
      int fp_e_12_13 = data_edge_firstparent_face[ e_12_13 ];
      int fp_e_12_23 = data_edge_firstparent_face[ e_12_23 ];
      int fp_e_13_23 = data_edge_firstparent_face[ e_13_23 ];
      
      int fp_e_02_13 = data_edge_firstparent_face[ e_02_13 ];
      
      
      
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      
      data_edge_firstparent_face[ e_01_02 ] = f_01_02_03;
      data_face_nextparents_of_edges[ f_01_02_03 ][ 0 ] = f_01_02_13; 
      data_face_nextparents_of_edges[ f_01_02_13 ][ 0 ] = fp_e_01_02; 
      
      data_edge_firstparent_face[ e_01_12 ] = f_01_12_13;
      data_face_nextparents_of_edges[ f_01_12_13 ][ 0 ] = fp_e_01_12;
      
      data_edge_firstparent_face[ e_02_12 ] = f_02_12_23;
      data_face_nextparents_of_edges[ f_02_12_23 ][ 0 ] = f_02_12_13;
      data_face_nextparents_of_edges[ f_02_12_13 ][ 0 ] = fp_e_02_12;
      
      
      data_edge_firstparent_face[ e_01_03 ] = f_01_02_03;
      data_face_nextparents_of_edges[ f_01_02_03 ][ 1 ] = fp_e_01_03;
      
      data_edge_firstparent_face[ e_01_13 ] = f_01_12_13;
      data_face_nextparents_of_edges[ f_01_12_13 ][ 1 ] = f_01_02_13;
      data_face_nextparents_of_edges[ f_01_02_13 ][ 1 ] = fp_e_01_13;
      
      data_edge_firstparent_face[ e_03_13 ] = f_03_13_23;
      data_face_nextparents_of_edges[ f_03_13_23 ][ 0 ] = f_02_03_13;
      data_face_nextparents_of_edges[ f_02_03_13 ][ 2 ] = fp_e_03_13;
      
      
      data_edge_firstparent_face[ e_02_03 ] = f_01_02_03;
      data_face_nextparents_of_edges[ f_01_02_03 ][ 2 ] = f_02_03_13;
      data_face_nextparents_of_edges[ f_02_03_13 ][ 0 ] = fp_e_02_03;
      
      data_edge_firstparent_face[ e_02_23 ] = f_02_12_23;
      data_face_nextparents_of_edges[ f_02_12_23 ][ 1 ] = f_02_13_23;
      data_face_nextparents_of_edges[ f_02_13_23 ][ 1 ] = fp_e_02_23;
      
      data_edge_firstparent_face[ e_03_23 ] = f_03_13_23;
      data_face_nextparents_of_edges[ f_03_13_23 ][ 1 ] = fp_e_03_23;
      
      
      data_edge_firstparent_face[ e_12_13 ] = f_01_12_13;
      data_face_nextparents_of_edges[ f_01_12_13 ][ 2 ] = f_02_12_13;
      data_face_nextparents_of_edges[ f_02_12_13 ][ 2 ] = fp_e_12_13;
      
      data_edge_firstparent_face[ e_12_23 ] = f_02_12_23;
      data_face_nextparents_of_edges[ f_02_12_23 ][ 2 ] = fp_e_12_23;
      
      data_edge_firstparent_face[ e_13_23 ] = f_03_13_23;
      data_face_nextparents_of_edges[ f_03_13_23 ][ 2 ] = f_02_13_23;
      data_face_nextparents_of_edges[ f_02_13_23 ][ 2 ] = fp_e_13_23;
      
      
      data_edge_firstparent_face[ e_02_13 ] = f_01_02_13;
      data_face_nextparents_of_edges[ f_01_02_13 ][ 2 ] = f_02_03_13;
      data_face_nextparents_of_edges[ f_02_03_13 ][ 1 ] = f_02_12_13;
      data_face_nextparents_of_edges[ f_02_12_13 ][ 1 ] = f_02_13_23;
      data_face_nextparents_of_edges[ f_02_13_23 ][ 0 ] = fp_e_02_13;
      
    }
    
    
    
    
    
    
    
    /*** SET THE VERTICES AND EDGES OF ALL FACES ****/
    
    /* for each new outer face, set the new vertices */
    
    for( int f = 0; f < counter_faces; f++ )
    {
      int v00 = data_face_vertices[f][0];
      int v11 = data_face_vertices[f][1];
      int v22 = data_face_vertices[f][2];
      
      int v01 = counter_vertices + data_face_edges[f][0];
      int v02 = counter_vertices + data_face_edges[f][1];
      int v12 = counter_vertices + data_face_edges[f][2];
      
      int f_00_01_12 = 0 * counter_faces + f;
      int f_01_11_12 = 1 * counter_faces + f;
      int f_02_12_22 = 2 * counter_faces + f;
      int f_01_02_12 = 3 * counter_faces + f;
      
        
      
      
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
      data_face_vertices[ f_00_01_12 ][0] = v00;
      data_face_vertices[ f_00_01_12 ][1] = v01;
      data_face_vertices[ f_00_01_12 ][2] = v02;
      
      data_face_vertices[ f_01_11_12 ][0] = v01;
      data_face_vertices[ f_01_11_12 ][1] = v11;
      data_face_vertices[ f_01_11_12 ][2] = v12;
      
      data_face_vertices[ f_02_12_22 ][0] = v02;
      data_face_vertices[ f_02_12_22 ][1] = v12;
      data_face_vertices[ f_02_12_22 ][2] = v22;
      
      data_face_vertices[ f_01_02_12 ][0] = v01;
      data_face_vertices[ f_01_02_12 ][1] = v02;
      data_face_vertices[ f_01_02_12 ][2] = v12;
      
    }
    
    /* for each new interior face, set the new vertices */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
      int v01 = counter_vertices + data_tetrahedron_edges[t][0];
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v03 = counter_vertices + data_tetrahedron_edges[t][2];
      int v12 = counter_vertices + data_tetrahedron_edges[t][3];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      int v23 = counter_vertices + data_tetrahedron_edges[t][5];
      
      assert( counter_vertices <= v01 && v01 < counter_vertices + counter_edges );
      assert( counter_vertices <= v02 && v02 < counter_vertices + counter_edges );
      assert( counter_vertices <= v03 && v03 < counter_vertices + counter_edges );
      assert( counter_vertices <= v12 && v12 < counter_vertices + counter_edges );
      assert( counter_vertices <= v13 && v13 < counter_vertices + counter_edges );
      assert( counter_vertices <= v23 && v23 < counter_vertices + counter_edges );
      
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      data_face_vertices[ f_01_02_03 ][0] = v01;
      data_face_vertices[ f_01_02_03 ][1] = v02;
      data_face_vertices[ f_01_02_03 ][2] = v03;
      
      data_face_vertices[ f_01_12_13 ][0] = v01;
      data_face_vertices[ f_01_12_13 ][1] = v12;
      data_face_vertices[ f_01_12_13 ][2] = v13;
      
      data_face_vertices[ f_02_12_23 ][0] = v02;
      data_face_vertices[ f_02_12_23 ][1] = v12;
      data_face_vertices[ f_02_12_23 ][2] = v23;
      
      data_face_vertices[ f_03_13_23 ][0] = v03;
      data_face_vertices[ f_03_13_23 ][1] = v13;
      data_face_vertices[ f_03_13_23 ][2] = v23;
      
      data_face_vertices[ f_01_02_13 ][0] = v01;
      data_face_vertices[ f_01_02_13 ][1] = v02;
      data_face_vertices[ f_01_02_13 ][2] = v13;
      
      data_face_vertices[ f_02_03_13 ][0] = v02;
      data_face_vertices[ f_02_03_13 ][1] = v03;
      data_face_vertices[ f_02_03_13 ][2] = v13;
      
      data_face_vertices[ f_02_12_13 ][0] = v02;
      data_face_vertices[ f_02_12_13 ][1] = v12;
      data_face_vertices[ f_02_12_13 ][2] = v13;
      
      data_face_vertices[ f_02_13_23 ][0] = v02;
      data_face_vertices[ f_02_13_23 ][1] = v13;
      data_face_vertices[ f_02_13_23 ][2] = v23;
      
      
      
    }
    
    /* for each new outer face, set the new edges */
    
    for( int f = 0; f < counter_faces; f++ )
    {
    
      int e00_01 = 0 * counter_edges + data_face_edges[f][0];
      int e00_02 = 0 * counter_edges + data_face_edges[f][1];
      int e01_11 = 1 * counter_edges + data_face_edges[f][0];
      int e02_22 = 1 * counter_edges + data_face_edges[f][1];
      int e11_12 = 0 * counter_edges + data_face_edges[f][2];
      int e12_22 = 1 * counter_edges + data_face_edges[f][2];
      
      int e01_02 = 2 * counter_edges + 0 * counter_faces + f;
      int e01_12 = 2 * counter_edges + 1 * counter_faces + f;
      int e02_12 = 2 * counter_edges + 2 * counter_faces + f;
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
      data_face_edges[ 0 * counter_faces + f ][0] = e00_01;
      data_face_edges[ 0 * counter_faces + f ][1] = e00_02;
      data_face_edges[ 0 * counter_faces + f ][2] = e01_02;
      
      data_face_edges[ 1 * counter_faces + f ][0] = e01_11;
      data_face_edges[ 1 * counter_faces + f ][1] = e01_12;
      data_face_edges[ 1 * counter_faces + f ][2] = e11_12;
      
      data_face_edges[ 2 * counter_faces + f ][0] = e02_12;
      data_face_edges[ 2 * counter_faces + f ][1] = e02_22;
      data_face_edges[ 2 * counter_faces + f ][2] = e12_22;
      
      data_face_edges[ 3 * counter_faces + f ][0] = e01_02;
      data_face_edges[ 3 * counter_faces + f ][1] = e01_12;
      data_face_edges[ 3 * counter_faces + f ][2] = e02_12;
      
    }
    
    
    /* for each new interior face, set the new edges */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      int f_012 = data_tetrahedron_faces[t][0];
      int f_013 = data_tetrahedron_faces[t][1];
      int f_023 = data_tetrahedron_faces[t][2];
      int f_123 = data_tetrahedron_faces[t][3];
      
      int e_01_02 = 2 * counter_edges + 0 * counter_faces + f_012;
      int e_01_12 = 2 * counter_edges + 1 * counter_faces + f_012;
      int e_02_12 = 2 * counter_edges + 2 * counter_faces + f_012;
      
      int e_01_03 = 2 * counter_edges + 0 * counter_faces + f_013;
      int e_01_13 = 2 * counter_edges + 1 * counter_faces + f_013;
      int e_03_13 = 2 * counter_edges + 2 * counter_faces + f_013;
      
      int e_02_03 = 2 * counter_edges + 0 * counter_faces + f_023;
      int e_02_23 = 2 * counter_edges + 1 * counter_faces + f_023;
      int e_03_23 = 2 * counter_edges + 2 * counter_faces + f_023;
      
      int e_12_13 = 2 * counter_edges + 0 * counter_faces + f_123;
      int e_12_23 = 2 * counter_edges + 1 * counter_faces + f_123;
      int e_13_23 = 2 * counter_edges + 2 * counter_faces + f_123;
      
      int e_02_13 = 2 * counter_edges + 3 * counter_faces + t;
      
      data_face_edges[ f_01_02_03 ][0] = e_01_02;
      data_face_edges[ f_01_02_03 ][1] = e_01_03;
      data_face_edges[ f_01_02_03 ][2] = e_02_03;
      
      data_face_edges[ f_01_12_13 ][0] = e_01_12;
      data_face_edges[ f_01_12_13 ][1] = e_01_13;
      data_face_edges[ f_01_12_13 ][2] = e_12_13;
      
      data_face_edges[ f_02_12_23 ][0] = e_02_12;
      data_face_edges[ f_02_12_23 ][1] = e_02_23;
      data_face_edges[ f_02_12_23 ][2] = e_12_23;
      
      data_face_edges[ f_03_13_23 ][0] = e_03_13;
      data_face_edges[ f_03_13_23 ][1] = e_03_23;
      data_face_edges[ f_03_13_23 ][2] = e_13_23;
      
      data_face_edges[ f_01_02_13 ][0] = e_01_02;
      data_face_edges[ f_01_02_13 ][1] = e_01_13;
      data_face_edges[ f_01_02_13 ][2] = e_02_13;
      
      data_face_edges[ f_02_03_13 ][0] = e_02_03;
      data_face_edges[ f_02_03_13 ][1] = e_02_13;
      data_face_edges[ f_02_03_13 ][2] = e_03_13;
      
      data_face_edges[ f_02_12_13 ][0] = e_02_12;
      data_face_edges[ f_02_12_13 ][1] = e_02_13;
      data_face_edges[ f_02_12_13 ][2] = e_12_13;
      
      data_face_edges[ f_02_13_23 ][0] = e_02_13;
      data_face_edges[ f_02_13_23 ][1] = e_02_23;
      data_face_edges[ f_02_13_23 ][2] = e_13_23;
      
    }
    
    
    
    /**************************/
    /***   ALL DIMENSIONS  ****/
    /**************************/
    
    
    /*** TREAT THE OLD VERTICES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    /* 1. for each old vertex, set the new parent tetrahedron */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      
      int p = data_vertex_firstparent_tetrahedron[v];
      
      assert( p != nullindex && 0 <= p && p < counter_tetrahedra );
      
      int vi = nullindex;
      if( data_tetrahedron_vertices[p][0] == v ) vi = 0;
      if( data_tetrahedron_vertices[p][1] == v ) vi = 1;
      if( data_tetrahedron_vertices[p][2] == v ) vi = 2;
      if( data_tetrahedron_vertices[p][3] == v ) vi = 3;
      assert( vi != nullindex );
      
      assert( data_tetrahedron_vertices[p][vi] == v );
      
      data_vertex_firstparent_tetrahedron[v] = vi * counter_tetrahedra + p;
      
    }
    
    
    /* 2. for each old tetrahedron, relocate the data of the old vertices' parent tetrahedron */
    
    for( int t  = 0; t  < counter_tetrahedra;  t++ )
    for( int vi = 0; vi <                  4; vi++ )
    {
      
      int q = data_tetrahedron_nextparents_of_vertices[t][vi];
      
      int v = data_tetrahedron_vertices[t][vi];
      
      if( q == nullindex ) {
        
        data_tetrahedron_nextparents_of_vertices[ vi * counter_tetrahedra + t ][vi] = nullindex;
        
      } else {
        
        int vinp = nullindex;
        if( data_tetrahedron_vertices[q][0] == v ) vinp = 0;
        if( data_tetrahedron_vertices[q][1] == v ) vinp = 1;
        if( data_tetrahedron_vertices[q][2] == v ) vinp = 2;
        if( data_tetrahedron_vertices[q][3] == v ) vinp = 3;
        assert( vinp != nullindex );
        
        assert( data_tetrahedron_vertices[q][vinp] == v );
        
        data_tetrahedron_nextparents_of_vertices[ vi * counter_tetrahedra + t ][vi] = vinp * counter_tetrahedra + q;
      
      } 
      
    }
    
    /*
    *       00 01 02 03
    *       01 11 12 13
    *       02 12 22 23
    *       03 13 23 33
    *       
    *       01 02 03 13
    *       01 02 12 13 
    *       02 03 13 23
    *       02 12 13 23
    */
       
        
    /*** TREAT THE NEW VERTICES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        
        int v01 = counter_vertices + data_tetrahedron_edges[ t ][ 0 ];
        int v02 = counter_vertices + data_tetrahedron_edges[ t ][ 1 ];
        int v03 = counter_vertices + data_tetrahedron_edges[ t ][ 2 ];
        int v12 = counter_vertices + data_tetrahedron_edges[ t ][ 3 ];
        int v13 = counter_vertices + data_tetrahedron_edges[ t ][ 4 ];
        int v23 = counter_vertices + data_tetrahedron_edges[ t ][ 5 ];
        
        int fp_v01 = data_vertex_firstparent_tetrahedron[ v01 ];
        int fp_v02 = data_vertex_firstparent_tetrahedron[ v02 ];
        int fp_v03 = data_vertex_firstparent_tetrahedron[ v03 ];
        int fp_v12 = data_vertex_firstparent_tetrahedron[ v12 ];
        int fp_v13 = data_vertex_firstparent_tetrahedron[ v13 ];
        int fp_v23 = data_vertex_firstparent_tetrahedron[ v23 ];
        
        int t_00_01_02_03 = 0 * counter_tetrahedra + t;
        int t_01_11_12_13 = 1 * counter_tetrahedra + t;
        int t_02_12_22_23 = 2 * counter_tetrahedra + t;
        int t_03_13_23_33 = 3 * counter_tetrahedra + t;
        int t_01_02_03_13 = 4 * counter_tetrahedra + t;
        int t_01_02_12_13 = 5 * counter_tetrahedra + t;
        int t_02_03_13_23 = 6 * counter_tetrahedra + t;
        int t_02_12_13_23 = 7 * counter_tetrahedra + t;
        
        data_vertex_firstparent_tetrahedron[ v01 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_vertices[ t_00_01_02_03 ][ 1 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_11_12_13 ][ 0 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_03_13 ][ 0 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_12_13 ][ 0 ] = fp_v01;
        
        data_vertex_firstparent_tetrahedron[ v02 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_vertices[ t_00_01_02_03 ][ 2 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_22_23 ][ 0 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_03_13 ][ 1 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_12_13 ][ 1 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_03_13_23 ][ 0 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_13_23 ][ 0 ] = fp_v02;
        
        data_vertex_firstparent_tetrahedron[ v03 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_vertices[ t_00_01_02_03 ][ 3 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_vertices[ t_03_13_23_33 ][ 0 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_03_13 ][ 2 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_03_13_23 ][ 1 ] = fp_v03;
        
        data_vertex_firstparent_tetrahedron[ v12 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_11_12_13 ][ 2 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_22_23 ][ 1 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_12_13 ][ 2 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_13_23 ][ 1 ] = fp_v12;
        
        data_vertex_firstparent_tetrahedron[ v13 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_11_12_13 ][ 3 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_vertices[ t_03_13_23_33 ][ 1 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_03_13 ][ 3 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_12_13 ][ 3 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_03_13_23 ][ 2 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_13_23 ][ 2 ] = fp_v13;
        
        data_vertex_firstparent_tetrahedron[ v23 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_22_23 ][ 3 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_vertices[ t_03_13_23_33 ][ 2 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_03_13_23 ][ 3 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_13_23 ][ 3 ] = fp_v23;
        
    }
    
    /*** TREAT THE BISECTED EDGES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        
        int e_00_01 = data_tetrahedron_edges[ t ][ 0 ];
        int e_00_02 = data_tetrahedron_edges[ t ][ 1 ];
        int e_00_03 = data_tetrahedron_edges[ t ][ 2 ];
        int e_11_12 = data_tetrahedron_edges[ t ][ 3 ];
        int e_11_13 = data_tetrahedron_edges[ t ][ 4 ];
        int e_22_23 = data_tetrahedron_edges[ t ][ 5 ];
        
        int e_01_11 = counter_edges + data_tetrahedron_edges[ t ][ 0 ];
        int e_02_22 = counter_edges + data_tetrahedron_edges[ t ][ 1 ];
        int e_03_33 = counter_edges + data_tetrahedron_edges[ t ][ 2 ];
        int e_12_22 = counter_edges + data_tetrahedron_edges[ t ][ 3 ];
        int e_13_33 = counter_edges + data_tetrahedron_edges[ t ][ 4 ];
        int e_23_33 = counter_edges + data_tetrahedron_edges[ t ][ 5 ];
        
        
        int fp_e_00_01 = data_edge_firstparent_tetrahedron[ e_00_01 ];
        int fp_e_00_02 = data_edge_firstparent_tetrahedron[ e_00_02 ];
        int fp_e_00_03 = data_edge_firstparent_tetrahedron[ e_00_03 ];
        int fp_e_11_12 = data_edge_firstparent_tetrahedron[ e_11_12 ];
        int fp_e_11_13 = data_edge_firstparent_tetrahedron[ e_11_13 ];
        int fp_e_22_23 = data_edge_firstparent_tetrahedron[ e_22_23 ];
        
        int fp_e_01_11 = data_edge_firstparent_tetrahedron[ e_01_11 ];
        int fp_e_02_22 = data_edge_firstparent_tetrahedron[ e_02_22 ];
        int fp_e_03_33 = data_edge_firstparent_tetrahedron[ e_03_33 ];
        int fp_e_12_22 = data_edge_firstparent_tetrahedron[ e_12_22 ];
        int fp_e_13_33 = data_edge_firstparent_tetrahedron[ e_13_33 ];
        int fp_e_23_33 = data_edge_firstparent_tetrahedron[ e_23_33 ];
        
        
        int t_00_01_02_03 = 0 * counter_tetrahedra + t;
        int t_01_11_12_13 = 1 * counter_tetrahedra + t;
        int t_02_12_22_23 = 2 * counter_tetrahedra + t;
        int t_03_13_23_33 = 3 * counter_tetrahedra + t;

        data_edge_firstparent_tetrahedron[ e_00_01 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 0 ] = fp_e_00_01;
        
        data_edge_firstparent_tetrahedron[ e_00_02 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 1 ] = fp_e_00_02;
        
        data_edge_firstparent_tetrahedron[ e_00_03 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 2 ] = fp_e_00_03;
        
        data_edge_firstparent_tetrahedron[ e_11_12 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 3 ] = fp_e_11_12;
        
        data_edge_firstparent_tetrahedron[ e_11_13 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 4 ] = fp_e_11_13;
        
        data_edge_firstparent_tetrahedron[ e_22_23 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 5 ] = fp_e_22_23;
        
        
        data_edge_firstparent_tetrahedron[ e_01_11 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 0 ] = fp_e_01_11;
        
        data_edge_firstparent_tetrahedron[ e_02_22 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 1 ] = fp_e_02_22;
        
        data_edge_firstparent_tetrahedron[ e_03_33 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 2 ] = fp_e_03_33;
        
        data_edge_firstparent_tetrahedron[ e_12_22 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 3 ] = fp_e_12_22;
        
        data_edge_firstparent_tetrahedron[ e_13_33 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 4 ] = fp_e_13_33;
        
        data_edge_firstparent_tetrahedron[ e_23_33 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 5 ] = fp_e_23_33;
                
        
    }
    
    
    
    
    /*** TREAT THE FACE-BASED EDGES AND THEIR CONNECTION TO TETRAHEDRA ****/
    /*** TREAT THE SINGLE INTERIOR EDGES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        
        int f012 = data_tetrahedron_faces[ t ][ 0 ];
        int f013 = data_tetrahedron_faces[ t ][ 1 ];
        int f023 = data_tetrahedron_faces[ t ][ 2 ];
        int f123 = data_tetrahedron_faces[ t ][ 2 ];
        
        /* 01 02 12 */
        /* 01 03 13 */
        /* 02 03 23 */
        /* 12 13 23 */
        
        int e_01_02 = 2 * counter_edges + 0 * counter_faces + f012;
        int e_01_12 = 2 * counter_edges + 1 * counter_faces + f012;
        int e_02_12 = 2 * counter_edges + 2 * counter_faces + f012;
        
        int e_01_03 = 2 * counter_edges + 0 * counter_faces + f013;
        int e_01_13 = 2 * counter_edges + 1 * counter_faces + f013;
        int e_03_13 = 2 * counter_edges + 2 * counter_faces + f013;
        
        int e_02_03 = 2 * counter_edges + 0 * counter_faces + f023;
        int e_02_23 = 2 * counter_edges + 1 * counter_faces + f023;
        int e_03_23 = 2 * counter_edges + 2 * counter_faces + f023;
        
        int e_12_13 = 2 * counter_edges + 0 * counter_faces + f123;
        int e_12_23 = 2 * counter_edges + 1 * counter_faces + f123;
        int e_13_23 = 2 * counter_edges + 2 * counter_faces + f123;
        
        int e_02_13 = 2 * counter_edges + 3 * counter_faces + t;
    
        int fp_e_01_02 = data_edge_firstparent_tetrahedron[ e_01_02 ];
        int fp_e_01_12 = data_edge_firstparent_tetrahedron[ e_01_12 ];
        int fp_e_02_12 = data_edge_firstparent_tetrahedron[ e_02_12 ];
        
        int fp_e_01_03 = data_edge_firstparent_tetrahedron[ e_01_03 ];
        int fp_e_01_13 = data_edge_firstparent_tetrahedron[ e_01_13 ];
        int fp_e_03_13 = data_edge_firstparent_tetrahedron[ e_03_13 ];
        
        int fp_e_02_03 = data_edge_firstparent_tetrahedron[ e_02_03 ];
        int fp_e_02_23 = data_edge_firstparent_tetrahedron[ e_02_23 ];
        int fp_e_03_23 = data_edge_firstparent_tetrahedron[ e_03_23 ];
        
        int fp_e_12_13 = data_edge_firstparent_tetrahedron[ e_12_13 ];
        int fp_e_12_23 = data_edge_firstparent_tetrahedron[ e_12_23 ];
        int fp_e_13_23 = data_edge_firstparent_tetrahedron[ e_13_23 ];
        
        int fp_e_02_13 = data_edge_firstparent_tetrahedron[ e_02_13 ];
        assert( fp_e_02_13 == nullindex );
        
        
        int t_00_01_02_03 = 0 * counter_tetrahedra + t;
        int t_01_11_12_13 = 1 * counter_tetrahedra + t;
        int t_02_12_22_23 = 2 * counter_tetrahedra + t;
        int t_03_13_23_33 = 3 * counter_tetrahedra + t;
        int t_01_02_03_13 = 4 * counter_tetrahedra + t;
        int t_01_02_12_13 = 5 * counter_tetrahedra + t;
        int t_02_03_13_23 = 6 * counter_tetrahedra + t;
        int t_02_12_13_23 = 7 * counter_tetrahedra + t;
        
        
        data_edge_firstparent_tetrahedron[ e_01_02 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 3 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_03_13 ][ 0 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 0 ] = fp_e_01_02;
        
        data_edge_firstparent_tetrahedron[ e_01_12 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 2 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 2 ] = fp_e_01_12;
        
        data_edge_firstparent_tetrahedron[ e_02_12 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 0 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 3 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 0 ] = fp_e_02_12;
        
        
        data_edge_firstparent_tetrahedron[ e_01_03 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 4 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_03_13 ][ 1 ] = fp_e_01_03;
        
        data_edge_firstparent_tetrahedron[ e_01_13 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 2 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_03_13 ][ 2 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 2 ] = fp_e_01_13;
        
        data_edge_firstparent_tetrahedron[ e_03_13 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 0 ] = fp_e_03_13;
        
        
        data_edge_firstparent_tetrahedron[ e_02_03 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 5 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_03_13 ][ 3 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_03_13_23 ][ 0 ] = fp_e_02_03;
        
        data_edge_firstparent_tetrahedron[ e_02_23 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 2 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_03_13_23 ][ 2 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 2 ] = fp_e_02_23;
        
        data_edge_firstparent_tetrahedron[ e_03_23 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 2 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_03_13_23 ][ 3 ] = fp_e_03_23;
        
        
        data_edge_firstparent_tetrahedron[ e_12_13 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 5 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 5 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 3 ] = fp_e_12_13;
        
        data_edge_firstparent_tetrahedron[ e_12_23 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 4 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 4 ] = fp_e_12_23;
        
        data_edge_firstparent_tetrahedron[ e_13_23 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 3 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_03_13_23 ][ 5 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 5 ] = fp_e_13_23;
        
        
        
        data_edge_firstparent_tetrahedron[ e_02_13 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_03_13 ][ 4 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 4 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_03_13_23 ][ 3 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 1 ] = fp_e_02_13;
        
        
        
        
    }
    
    
    
    
    
    /*** TREAT THE OUTER FACES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    /* TODO */
    
    
    /*** TREAT THE INNER FACES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    
    
    
    
    
    
    /* for each new tetrahedron, set the new vertices */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
      int v00 = data_tetrahedron_vertices[t][0];
      int v11 = data_tetrahedron_vertices[t][1];
      int v22 = data_tetrahedron_vertices[t][2];
      int v33 = data_tetrahedron_vertices[t][3];
      
      int v01 = counter_vertices + data_tetrahedron_edges[t][0];
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v03 = counter_vertices + data_tetrahedron_edges[t][2];
      int v12 = counter_vertices + data_tetrahedron_edges[t][3];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      int v23 = counter_vertices + data_tetrahedron_edges[t][5];
      
      //       00 01 02 03
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][0] = v00;
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][1] = v01;
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][2] = v02;
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][3] = v03;
      
      //       01 11 12 13
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][0] = v01;
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][1] = v11;
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][2] = v12;
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][3] = v13;
      
      //       02 12 22 23
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][0] = v02;
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][1] = v12;
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][2] = v22;
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][3] = v23;
      
      //       03 13 23 33
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][0] = v03;
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][1] = v13;
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][2] = v23;
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][3] = v33;
      
      
      //       01 02 03 13
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][0] = v01;
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][1] = v02;
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][2] = v03;
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][3] = v13;
      
      //       01 02 12 13 
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][0] = v01;
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][1] = v02;
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][2] = v12;
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][3] = v13;
      
      //       02 03 13 23
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][0] = v02;
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][1] = v03;
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][2] = v13;
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][3] = v23;
      
      //       02 12 13 23
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][0] = v02;
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][1] = v12;
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][2] = v13;
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][3] = v23;
      
    }
    
    /* for each new tetrahedron, set the new edges */
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        // 00 11 22 33 
        
        // From old edges ...
        int e00_01 = data_tetrahedron_edges[t][0] + 0 * counter_edges;
        int e01_11 = data_tetrahedron_edges[t][0] + 1 * counter_edges;
        
        int e00_02 = data_tetrahedron_edges[t][1] + 0 * counter_edges;
        int e02_22 = data_tetrahedron_edges[t][1] + 1 * counter_edges;
        
        int e00_03 = data_tetrahedron_edges[t][2] + 0 * counter_edges;
        int e03_33 = data_tetrahedron_edges[t][2] + 1 * counter_edges;
        
        int e11_12 = data_tetrahedron_edges[t][3] + 0 * counter_edges;
        int e12_22 = data_tetrahedron_edges[t][3] + 1 * counter_edges;
        
        int e11_13 = data_tetrahedron_edges[t][4] + 0 * counter_edges;
        int e13_33 = data_tetrahedron_edges[t][4] + 1 * counter_edges;
        
        int e22_23 = data_tetrahedron_edges[t][5] + 0 * counter_edges;
        int e23_33 = data_tetrahedron_edges[t][5] + 1 * counter_edges;
        
        
        // ... from old faces ... 
        
        // 00 11 22  
        int e01_02 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][0];
        int e01_12 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][0];
        int e02_12 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][0];
        
        // 00 11 33 
        int e01_03 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][1];
        int e01_13 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][1];
        int e03_13 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][1];
        
        // 00 22 33 
        int e02_03 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][2];
        int e02_23 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][2];
        int e03_23 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][2];
        
        // 11 22 33 
        int e12_13 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][3];
        int e12_23 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][3];
        int e13_23 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][3];
        
        // ... and the single new internal one. 
        
        int e02_13 = 2 * counter_edges + 3 * counter_faces + t;
        
        
        // fill in to the tetrahedra
        
        //       00 01 02 03
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][0] = e00_01;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][1] = e00_02;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][2] = e00_03;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][3] = e01_02;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][4] = e01_03;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][5] = e02_03;
        
        //       01 11 12 13
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][0] = e01_11;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][1] = e01_12;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][2] = e01_13;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][3] = e11_12;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][4] = e11_13;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][5] = e12_13;
        
        //       02 12 22 23
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][0] = e02_12;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][1] = e02_22;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][2] = e02_23;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][3] = e12_22;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][4] = e12_23;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][5] = e22_23;
        
        //       03 13 23 33
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][0] = e03_13;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][1] = e03_23;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][2] = e03_33;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][3] = e13_23;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][4] = e13_33;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][5] = e23_33;
        
        
        //       01 02 03 13
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][0] = e01_02;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][1] = e01_03;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][2] = e01_13;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][3] = e02_03;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][4] = e02_13;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][5] = e03_13;
        
        //       01 02 12 13 
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][0] = e01_02;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][1] = e01_12;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][2] = e01_13;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][3] = e02_12;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][4] = e02_13;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][5] = e12_13;
        
        //       02 03 13 23
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][0] = e02_03;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][1] = e02_13;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][2] = e02_23;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][3] = e03_13;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][4] = e03_23;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][5] = e13_23;
        
        //       02 12 13 23
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][0] = e02_12;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][1] = e02_13;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][2] = e02_23;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][3] = e12_13;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][4] = e12_23;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][5] = e13_23;
        
    }
    
    /* for each new tetrahedron, set the new faces */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      // 00 11 22 33 
      
      // From old faces ... 
      
      // 00 11 22 
      int f00_01_02 = data_tetrahedron_faces[t][0] + 0 * counter_faces;
      int f01_11_12 = data_tetrahedron_faces[t][0] + 1 * counter_faces;
      int f02_12_22 = data_tetrahedron_faces[t][0] + 2 * counter_faces;
      int f01_02_12 = data_tetrahedron_faces[t][0] + 3 * counter_faces;
      
      // 00 11 33 
      int f00_01_03 = data_tetrahedron_faces[t][1] + 0 * counter_faces;
      int f01_11_13 = data_tetrahedron_faces[t][1] + 1 * counter_faces;
      int f03_13_33 = data_tetrahedron_faces[t][1] + 2 * counter_faces;
      int f01_03_13 = data_tetrahedron_faces[t][1] + 3 * counter_faces;
      
      // 00 22 33 
      int f00_02_03 = data_tetrahedron_faces[t][2] + 0 * counter_faces;
      int f02_22_23 = data_tetrahedron_faces[t][2] + 1 * counter_faces;
      int f03_23_33 = data_tetrahedron_faces[t][2] + 2 * counter_faces;
      int f02_03_23 = data_tetrahedron_faces[t][2] + 3 * counter_faces;
      
      // 11 22 33 
      int f11_12_13 = data_tetrahedron_faces[t][3] + 0 * counter_faces;
      int f12_22_23 = data_tetrahedron_faces[t][3] + 1 * counter_faces;
      int f13_23_33 = data_tetrahedron_faces[t][3] + 2 * counter_faces;
      int f12_13_23 = data_tetrahedron_faces[t][3] + 3 * counter_faces;
      
      // ... new internal faces. 
      
      // 01 02 03
      // 01 12 13
      // 02 12 23
      // 03 13 23
      int f01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      
      // 01 02 13
      // 02 03 13
      // 02 12 13
      // 02 13 23
      int f01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      
      
      
      //       00 01 02 03
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][0] = f00_01_02;
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][1] = f00_01_03;
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][2] = f00_02_03;
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][3] = f01_02_03;
      
      //       01 11 12 13
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][0] = f01_11_12;
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][1] = f01_11_13;
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][2] = f01_12_13;
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][3] = f11_12_13;
      
      //       02 12 22 23
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][0] = f02_12_22;
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][1] = f02_12_23;
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][2] = f02_22_23;
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][3] = f12_22_23;
      
      //       03 13 23 33
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][0] = f03_13_23;
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][1] = f03_13_33;
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][2] = f03_23_33;
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][3] = f13_23_33;
      
      
      //       01 02 03 13
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][0] = f01_02_03;
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][1] = f01_02_13;
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][2] = f01_03_13;
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][3] = f02_03_13;
      
      //       01 02 12 13 
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][0] = f01_02_12;
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][1] = f01_02_13;
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][2] = f01_12_13;
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][3] = f02_12_13;
      
      //       02 03 13 23
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][0] = f02_03_13;
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][1] = f02_03_23;
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][2] = f02_13_23;
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][3] = f03_13_23;
      
      //       02 12 13 23
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][0] = f02_12_13;
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][1] = f02_12_23;
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][2] = f02_13_23;
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][3] = f12_13_23;
      
    }
    
    
    
    
    /* update the counters */
    
    counter_vertices   = new_counter_vertices;
    counter_edges      = new_counter_edges;
    counter_faces      = new_counter_faces;
    counter_tetrahedra = new_counter_tetrahedra;
    
    /* DONE */
    
//     check();
}













void MeshSimplicial3D::midpoint_refinement( int t )
{
    check();
    
    assert( 0 <= t && t < counter_faces );
    
    /* Allocate memory */
    
    data_tetrahedron_faces.resize               ( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
    data_face_firstparent_tetrahedron.resize    ( counter_faces + 6,                                         nullindex   );
    data_tetrahedron_nextparents_of_faces.resize( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
    
    data_tetrahedron_edges.resize               ( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    data_edge_firstparent_tetrahedron.resize    ( counter_edges + 4,                                                               nullindex   );
    data_tetrahedron_nextparents_of_edges.resize( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    
    data_tetrahedron_vertices.resize               ( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
    data_vertex_firstparent_tetrahedron.resize     ( counter_vertices + 1,                                      nullindex   );
    data_tetrahedron_nextparents_of_vertices.resize( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
    
    data_face_edges.resize               ( counter_faces + 6, { nullindex, nullindex, nullindex } );
    data_edge_firstparent_face.resize    ( counter_edges + 4,                         nullindex   );
    data_face_nextparents_of_edges.resize( counter_faces + 6, { nullindex, nullindex, nullindex } );
    
    data_face_vertices.resize               ( counter_faces + 6, { nullindex, nullindex, nullindex } );
    data_vertex_firstparent_face.resize     ( counter_vertices + 1,                      nullindex   );
    data_face_nextparents_of_vertices.resize( counter_faces + 6, { nullindex, nullindex, nullindex } );
    
    data_edge_vertices.resize               ( counter_edges + 4,   { nullindex, nullindex } );
    data_vertex_firstparent_edge.resize     ( counter_vertices + 1,             nullindex   );
    data_edge_nextparents_of_vertices.resize( counter_edges + 4,   { nullindex, nullindex } );
    
    getcoordinates().addcoordinates( 1 );
    
    
    /* load the new coordinate */
    
    getcoordinates().loadvector( counter_vertices, get_face_midpoint( t ) );
    
    
    /* assemble the data and auxiliary variables */
    int t_f0 = data_tetrahedron_faces[t][0];
    int t_f1 = data_tetrahedron_faces[t][1];
    int t_f2 = data_tetrahedron_faces[t][2];
    int t_f3 = data_tetrahedron_faces[t][3];
    
    int t_e0 = data_tetrahedron_edges[t][0];
    int t_e1 = data_tetrahedron_edges[t][1];
    int t_e2 = data_tetrahedron_edges[t][2];
    int t_e3 = data_tetrahedron_edges[t][3];
    int t_e4 = data_tetrahedron_edges[t][4];
    int t_e5 = data_tetrahedron_edges[t][5];
    
    int t_v0 = data_tetrahedron_vertices[t][0];
    int t_v1 = data_tetrahedron_vertices[t][1];
    int t_v2 = data_tetrahedron_vertices[t][2];
    int t_v3 = data_tetrahedron_vertices[t][3];
    
    int t_f_n0 = data_tetrahedron_nextparents_of_faces[t][0];
    int t_f_n1 = data_tetrahedron_nextparents_of_faces[t][1];
    int t_f_n2 = data_tetrahedron_nextparents_of_faces[t][2];
    int t_f_n3 = data_tetrahedron_nextparents_of_faces[t][3];
    
    int t_e_n0 = data_tetrahedron_nextparents_of_edges[t][0];
    int t_e_n1 = data_tetrahedron_nextparents_of_edges[t][1];
    int t_e_n2 = data_tetrahedron_nextparents_of_edges[t][2];
    int t_e_n3 = data_tetrahedron_nextparents_of_edges[t][3];
    int t_e_n4 = data_tetrahedron_nextparents_of_edges[t][4];
    int t_e_n5 = data_tetrahedron_nextparents_of_edges[t][5];
    
    int t_v_n0 = data_tetrahedron_nextparents_of_vertices[t][0];
    int t_v_n1 = data_tetrahedron_nextparents_of_vertices[t][1];
    int t_v_n2 = data_tetrahedron_nextparents_of_vertices[t][2];
    int t_v_n3 = data_tetrahedron_nextparents_of_vertices[t][3];
    
    
    int t0 = t;
    int t1 = counter_faces;
    int t2 = counter_faces + 1;
    int t3 = counter_faces + 2;
    
    int f0 = counter_faces + 0;
    int f1 = counter_faces + 1;
    int f2 = counter_faces + 2;
    int f3 = counter_faces + 3;
    int f4 = counter_faces + 4;
    int f5 = counter_faces + 5;
    
    int e0 = counter_edges + 0;
    int e1 = counter_edges + 1;
    int e2 = counter_edges + 2;
    int e3 = counter_edges + 3;
    
    int vn = counter_vertices;
    
    
    /* fill in: downward */ // TODO
    data_tetrahedron_vertices[ t0 ][0] = vn;
    data_tetrahedron_vertices[ t0 ][1] = t_v0;
    data_tetrahedron_vertices[ t0 ][2] = t_v1;
    data_tetrahedron_vertices[ t0 ][3] = t_v2;
    
    data_tetrahedron_vertices[ t1 ][0] = vn;
    data_tetrahedron_vertices[ t1 ][1] = t_v0;
    data_tetrahedron_vertices[ t1 ][2] = t_v1;
    data_tetrahedron_vertices[ t1 ][3] = t_v3;
    
    data_tetrahedron_vertices[ t2 ][0] = vn;
    data_tetrahedron_vertices[ t2 ][1] = t_v0;
    data_tetrahedron_vertices[ t2 ][2] = t_v2;
    data_tetrahedron_vertices[ t2 ][3] = t_v3;
    
    data_tetrahedron_vertices[ t3 ][0] = vn;
    data_tetrahedron_vertices[ t3 ][1] = t_v1;
    data_tetrahedron_vertices[ t3 ][2] = t_v2;
    data_tetrahedron_vertices[ t3 ][3] = t_v3;
    
    data_tetrahedron_edges[ t0 ][0] = e0;
    data_tetrahedron_edges[ t0 ][1] = e1;
    data_tetrahedron_edges[ t0 ][2] = e2;
    data_tetrahedron_edges[ t0 ][3] = t_e0;
    data_tetrahedron_edges[ t0 ][4] = t_e1;
    data_tetrahedron_edges[ t0 ][5] = t_e3;
    
    data_tetrahedron_edges[ t1 ][0] = e0;
    data_tetrahedron_edges[ t1 ][1] = e1;
    data_tetrahedron_edges[ t1 ][2] = e3;
    data_tetrahedron_edges[ t1 ][3] = t_e0;
    data_tetrahedron_edges[ t1 ][4] = t_e2;
    data_tetrahedron_edges[ t1 ][5] = t_e4;
    
    data_tetrahedron_edges[ t2 ][0] = e0;
    data_tetrahedron_edges[ t2 ][1] = e2;
    data_tetrahedron_edges[ t2 ][2] = e3;
    data_tetrahedron_edges[ t2 ][3] = t_e1;
    data_tetrahedron_edges[ t2 ][4] = t_e2;
    data_tetrahedron_edges[ t2 ][5] = t_e5;
    
    data_tetrahedron_edges[ t3 ][0] = e1;
    data_tetrahedron_edges[ t3 ][1] = e2;
    data_tetrahedron_edges[ t3 ][2] = e3;
    data_tetrahedron_edges[ t3 ][3] = t_e3;
    data_tetrahedron_edges[ t3 ][4] = t_e4;
    data_tetrahedron_edges[ t3 ][5] = t_e5;
    
    data_tetrahedron_faces[ t0 ][0] = f0;
    data_tetrahedron_faces[ t0 ][1] = f1;
    data_tetrahedron_faces[ t0 ][2] = f3;
    data_tetrahedron_faces[ t0 ][3] = t_f0;
    
    data_tetrahedron_faces[ t1 ][0] = f0;
    data_tetrahedron_faces[ t1 ][1] = f2;
    data_tetrahedron_faces[ t1 ][2] = f4;
    data_tetrahedron_faces[ t1 ][3] = t_f1;
    
    data_tetrahedron_faces[ t2 ][0] = f1;
    data_tetrahedron_faces[ t2 ][1] = f2;
    data_tetrahedron_faces[ t2 ][2] = f5;
    data_tetrahedron_faces[ t2 ][3] = t_f2;
    
    data_tetrahedron_faces[ t3 ][0] = f3;
    data_tetrahedron_faces[ t3 ][1] = f4;
    data_tetrahedron_faces[ t3 ][2] = f5;
    data_tetrahedron_faces[ t3 ][3] = t_f3;
    
    data_face_vertices[ f0 ][0] = vn;
    data_face_vertices[ f0 ][1] = t_v0;
    data_face_vertices[ f0 ][2] = t_v1;
    
    data_face_vertices[ f1 ][0] = vn;
    data_face_vertices[ f1 ][1] = t_v0;
    data_face_vertices[ f1 ][2] = t_v2;
    
    data_face_vertices[ f2 ][0] = vn;
    data_face_vertices[ f2 ][1] = t_v0;
    data_face_vertices[ f2 ][2] = t_v3;
    
    data_face_vertices[ f3 ][0] = vn;
    data_face_vertices[ f3 ][1] = t_v1;
    data_face_vertices[ f3 ][2] = t_v2;
    
    data_face_vertices[ f4 ][0] = vn;
    data_face_vertices[ f4 ][1] = t_v1;
    data_face_vertices[ f4 ][2] = t_v3;
    
    data_face_vertices[ f5 ][0] = vn;
    data_face_vertices[ f5 ][1] = t_v2;
    data_face_vertices[ f5 ][2] = t_v3;
    
    data_face_edges[ f0 ][0] = e0;
    data_face_edges[ f0 ][1] = e1;
    data_face_edges[ f0 ][2] = t_e0;
    
    data_face_edges[ f1 ][0] = e0;
    data_face_edges[ f1 ][1] = e2;
    data_face_edges[ f1 ][2] = t_e1;
    
    data_face_edges[ f2 ][0] = e0;
    data_face_edges[ f2 ][1] = e3;
    data_face_edges[ f2 ][2] = t_e2;
    
    data_face_edges[ f3 ][0] = e1;
    data_face_edges[ f3 ][1] = e2;
    data_face_edges[ f3 ][2] = t_e3;
    
    data_face_edges[ f4 ][0] = e1;
    data_face_edges[ f4 ][1] = e3;
    data_face_edges[ f4 ][2] = t_e4;
    
    data_face_edges[ f5 ][0] = e2;
    data_face_edges[ f5 ][1] = e3;
    data_face_edges[ f5 ][2] = t_e5;
    
    data_edge_vertices[ e0 ][0] = vn;
    data_edge_vertices[ e0 ][1] = t_v0;
    
    data_edge_vertices[ e1 ][0] = vn;
    data_edge_vertices[ e1 ][1] = t_v1;
    
    data_edge_vertices[ e2 ][0] = vn;
    data_edge_vertices[ e2 ][1] = t_v2;
    
    data_edge_vertices[ e3 ][0] = vn;
    data_edge_vertices[ e3 ][1] = t_v3;
    
    /* fill in: parentlist */
    data_tetrahedron_nextparents_of_vertices[ t0 ][0] = t1;
    data_tetrahedron_nextparents_of_vertices[ t0 ][1] = t1;
    data_tetrahedron_nextparents_of_vertices[ t0 ][2] = t1;
    data_tetrahedron_nextparents_of_vertices[ t0 ][3] = t2;
    
    data_tetrahedron_nextparents_of_vertices[ t1 ][0] = t2;
    data_tetrahedron_nextparents_of_vertices[ t1 ][1] = t2;
    data_tetrahedron_nextparents_of_vertices[ t1 ][2] = t3;
    data_tetrahedron_nextparents_of_vertices[ t1 ][3] = t2;
    
    data_tetrahedron_nextparents_of_vertices[ t2 ][0] = t3;
    data_tetrahedron_nextparents_of_vertices[ t2 ][1] = data_vertex_firstparent_tetrahedron[ t_v0 ];
    data_tetrahedron_nextparents_of_vertices[ t2 ][2] = t3;
    data_tetrahedron_nextparents_of_vertices[ t2 ][3] = t3;
    data_vertex_firstparent_tetrahedron[ t_v0 ] = t0;
    
    data_tetrahedron_nextparents_of_vertices[ t3 ][0] = nullindex;
    data_tetrahedron_nextparents_of_vertices[ t3 ][1] = data_vertex_firstparent_tetrahedron[ t_v1 ];
    data_tetrahedron_nextparents_of_vertices[ t3 ][2] = data_vertex_firstparent_tetrahedron[ t_v2 ];
    data_tetrahedron_nextparents_of_vertices[ t3 ][3] = data_vertex_firstparent_tetrahedron[ t_v3 ];
    data_vertex_firstparent_tetrahedron[ t_v1 ] = t0;
    data_vertex_firstparent_tetrahedron[ t_v2 ] = t0;
    data_vertex_firstparent_tetrahedron[ t_v3 ] = t1;
    
    data_tetrahedron_nextparents_of_edges[ t0 ][0] = t1; // e0
    data_tetrahedron_nextparents_of_edges[ t0 ][1] = t1; // e1
    data_tetrahedron_nextparents_of_edges[ t0 ][2] = t2; // e2
    data_tetrahedron_nextparents_of_edges[ t0 ][3] = t1; // t_e_n0
    data_tetrahedron_nextparents_of_edges[ t0 ][4] = t2; // t_e_n1
    data_tetrahedron_nextparents_of_edges[ t0 ][5] = t3 // t_e_n3
    
    data_tetrahedron_nextparents_of_edges[ t1 ][0] = t2; // e0
    data_tetrahedron_nextparents_of_edges[ t1 ][1] = t3; // e1
    data_tetrahedron_nextparents_of_edges[ t1 ][2] = t2; // e3
    data_tetrahedron_nextparents_of_edges[ t1 ][3] = t_e_n0; //FIXME: relink below 
    data_tetrahedron_nextparents_of_edges[ t1 ][4] = t2; // t_e_n2
    data_tetrahedron_nextparents_of_edges[ t1 ][5] = t3; // t_e_n4
    
    data_tetrahedron_nextparents_of_edges[ t2 ][0] = nullindex; data_edge_firstparent_tetrahedron[ e0 ] = t0;
    data_tetrahedron_nextparents_of_edges[ t2 ][1] = t3; // e2
    data_tetrahedron_nextparents_of_edges[ t2 ][2] = t3; // e3
    data_tetrahedron_nextparents_of_edges[ t2 ][3] = t_e_n1; //FIXME: relink below
    data_tetrahedron_nextparents_of_edges[ t2 ][4] = t_e_n2; //FIXME: relink below
    data_tetrahedron_nextparents_of_edges[ t2 ][5] = t3; // t_e_n5
    
    data_tetrahedron_nextparents_of_edges[ t3 ][0] = nullindex; data_edge_firstparent_tetrahedron[ e1 ] = t0;
    data_tetrahedron_nextparents_of_edges[ t3 ][1] = nullindex; data_edge_firstparent_tetrahedron[ e2 ] = t0;
    data_tetrahedron_nextparents_of_edges[ t3 ][2] = nullindex; data_edge_firstparent_tetrahedron[ e3 ] = t0;
    data_tetrahedron_nextparents_of_edges[ t3 ][3] = t_e_n3; //FIXME: relink below
    data_tetrahedron_nextparents_of_edges[ t3 ][4] = t_e_n4; //FIXME: relink below
    data_tetrahedron_nextparents_of_edges[ t3 ][5] = t_e_n5; //FIXME: relink below
    
    data_tetrahedron_nextparents_of_faces[ t0 ][0] = t1;
    data_tetrahedron_nextparents_of_faces[ t0 ][1] = t2;
    data_tetrahedron_nextparents_of_faces[ t0 ][2] = t3;
    data_tetrahedron_nextparents_of_faces[ t0 ][3] = t_f_n0;
    
    data_tetrahedron_nextparents_of_faces[ t1 ][0] = nullindex;
    data_tetrahedron_nextparents_of_faces[ t1 ][1] = t2;
    data_tetrahedron_nextparents_of_faces[ t1 ][2] = t3;
    data_tetrahedron_nextparents_of_faces[ t1 ][3] = t_f_n1;
    
    data_tetrahedron_nextparents_of_faces[ t2 ][0] = nullindex;
    data_tetrahedron_nextparents_of_faces[ t2 ][1] = nullindex;
    data_tetrahedron_nextparents_of_faces[ t2 ][2] = t3;
    data_tetrahedron_nextparents_of_faces[ t2 ][3] = t_f_n2;
    
    data_tetrahedron_nextparents_of_faces[ t3 ][0] = nullindex;
    data_tetrahedron_nextparents_of_faces[ t3 ][1] = nullindex;
    data_tetrahedron_nextparents_of_faces[ t3 ][2] = nullindex;
    data_tetrahedron_nextparents_of_faces[ t3 ][3] = t_f_n3;
    
    data_face_nextparents_of_vertices[ f0 ][0] = f1;
    data_face_nextparents_of_vertices[ f0 ][1] = f1; // t_v0
    data_face_nextparents_of_vertices[ f0 ][2] = f3; // t_v1
    
    data_face_nextparents_of_vertices[ f1 ][0] = f2;
    data_face_nextparents_of_vertices[ f1 ][1] = f2; // t_v0
    data_face_nextparents_of_vertices[ f1 ][2] = f3; // t_v2
    
    data_face_nextparents_of_vertices[ f2 ][0] = f3;
    data_face_nextparents_of_vertices[ f2 ][1] = data_vertex_firstparent_face[ t_v0 ]; // t_v0
    data_face_nextparents_of_vertices[ f2 ][2] = f4; // t_v3
    
    data_face_nextparents_of_vertices[ f3 ][0] = f4;
    data_face_nextparents_of_vertices[ f3 ][1] = f4; // t_v1
    data_face_nextparents_of_vertices[ f3 ][2] = f5; // t_v2
    
    data_face_nextparents_of_vertices[ f4 ][0] = f5;
    data_face_nextparents_of_vertices[ f4 ][1] = data_vertex_firstparent_face[ t_v1 ]; // t_v1
    data_face_nextparents_of_vertices[ f4 ][2] = f5; // t_v3
    
    data_face_nextparents_of_vertices[ f5 ][0] = nullindex; data_vertex_firstparent_face[ vn ] = f0;
    data_face_nextparents_of_vertices[ f5 ][1] = data_vertex_firstparent_face[ t_v2 ]; // t_v2
    data_face_nextparents_of_vertices[ f5 ][2] = data_vertex_firstparent_face[ t_v3 ]; // t_v3
    
    data_vertex_firstparent_face[ t_v0 ] = f0;
    data_vertex_firstparent_face[ t_v1 ] = f0;
    data_vertex_firstparent_face[ t_v2 ] = f1;
    data_vertex_firstparent_face[ t_v3 ] = f2;
    
    data_face_nextparents_of_edges[ f0 ][0] = f1; // e0
    data_face_nextparents_of_edges[ f0 ][1] = f3; // e1
    data_face_nextparents_of_edges[ f0 ][2] = data_edge_firstparent_face[ t_e0 ]; // t_e0
    data_edge_firstparent_face[ t_e0 ] = f0;
    
    data_face_nextparents_of_edges[ f1 ][0] = f2; // e0
    data_face_nextparents_of_edges[ f1 ][1] = f3; // e2
    data_face_nextparents_of_edges[ f1 ][2] = data_edge_firstparent_face[ t_e1 ]; // t_e1
    data_edge_firstparent_face[ t_e1 ] = f1;
    
    data_face_nextparents_of_edges[ f2 ][0] = nullindex; // e0
    data_face_nextparents_of_edges[ f2 ][1] = f4; // e3
    data_face_nextparents_of_edges[ f2 ][2] = data_edge_firstparent_face[ t_e2 ]; // t_e2
    data_edge_firstparent_face[ t_e2 ] = f2;
    
    data_face_nextparents_of_edges[ f3 ][0] = f4; // e1
    data_face_nextparents_of_edges[ f3 ][1] = f5; // e2
    data_face_nextparents_of_edges[ f3 ][2] = data_edge_firstparent_face[ t_e3 ]; // t_e3
    data_edge_firstparent_face[ t_e3 ] = f3;
    
    data_face_nextparents_of_edges[ f4 ][0] = nullindex; // e1
    data_face_nextparents_of_edges[ f4 ][1] = f5; // e3
    data_face_nextparents_of_edges[ f4 ][2] = data_edge_firstparent_face[ t_e4 ]; // t_e4
    data_edge_firstparent_face[ t_e4 ] = f4;
    
    data_face_nextparents_of_edges[ f5 ][0] = nullindex; // e2
    data_face_nextparents_of_edges[ f5 ][1] = nullindex; // e3
    data_face_nextparents_of_edges[ f5 ][2] = data_edge_firstparent_face[ t_e5 ]; // t_e5
    data_edge_firstparent_face[ t_e5 ] = f5;
    
    data_edge_firstparent_face[ e0 ] = f0;
    data_edge_firstparent_face[ e1 ] = f0;
    data_edge_firstparent_face[ e2 ] = f1;
    data_edge_firstparent_face[ e3 ] = f4;
    
    
    data_edge_nextparents_of_vertices[ e0 ][0] = t_v0;
    data_edge_nextparents_of_vertices[ e0 ][1] = vn;
    
    data_edge_nextparents_of_vertices[ e1 ][0] = t_v1;
    data_edge_nextparents_of_vertices[ e1 ][1] = vn;
    
    data_edge_nextparents_of_vertices[ e2 ][0] = vn;
    data_edge_nextparents_of_vertices[ e2 ][1] = t_v2;
    
    data_edge_nextparents_of_vertices[ e3 ][0] = vn;
    data_edge_nextparents_of_vertices[ e3 ][1] = t_v2;
    
    
    /* face t_f0: nothing needs to change */
    
    /* face t_f1: relink */
    if( data_face_firstparent_tetrahedron[ t_f1 ] == t ) 
      data_face_firstparent_tetrahedron[ t_f1 ] = t1;
    else {
      int current_tetrahedron = data_face_firstparent_tetrahedron[ t_f1 ];
      while( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ] != t 
             &&
             data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ] != nullindex )
        current_tetrahedron = data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ];
      assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ] != nullindex );
      assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ] == t );
      data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ] = t1;
    }
    
    /* face t_f2: relink */
    if( data_face_firstparent_tetrahedron[ t_f2 ] == t ) 
      data_face_firstparent_tetrahedron[ t_f2 ] = t2;
    else {
      int current_tetrahedron = data_face_firstparent_tetrahedron[ t_f2 ];
      while( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ] != t 
             &&
             data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ] != nullindex )
        current_tetrahedron = data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ];
      assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ] != nullindex );
      assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ] == t );
      data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ] = t2;
    }
    
    /* face t_f3: relink */
    if( data_face_firstparent_tetrahedron[ t_f3 ] == t ) 
      data_face_firstparent_tetrahedron[ t_f3 ] = t3;
    else {
      int current_tetrahedron = data_face_firstparent_tetrahedron[ t_f3 ];
      while( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ] != t 
             &&
             data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ] != nullindex )
        current_tetrahedron = data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ];
      assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ] != nullindex );
      assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ] == t );
      data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ] = t3;
    }
    
    
    /* edge t_e0: nothing needs to change */
    
    /* edge t_e1: nothing needs to change */
    
    /* edge t_e2: relink */
    if( data_edge_firstparent_face[ t_e2 ] == t ) 
      data_edge_firstparent_face[ t_e2 ] = t1; // FIXME
    else {
      int current_face = data_edge_firstparent_face[ t_e2 ];
      while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ] != t 
             &&
             data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ] != nullindex )
        current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ];
      assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ] != nullindex );
      assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ] == t );
      data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ] = t1; // FIXME
    }
        
    /* edge t_e3: nothing needs to change */
    
    /* edge t_e4: relink */
    if( data_edge_firstparent_face[ t_e4 ] == t ) 
      data_edge_firstparent_face[ t_e4 ] = t4; // FIXME
    else {
      int current_face = data_edge_firstparent_face[ t_e4 ];
      while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ] != t 
             &&
             data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ] != nullindex )
        current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ];
      assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ] != nullindex );
      assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ] == t );
      data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ] = t4; // FIXME
    }
    
    /* edge t_e5: relink */
    if( data_edge_firstparent_face[ t_e5 ] == t ) 
      data_edge_firstparent_face[ t_e5 ] = t5; // FIXME
    else {
      int current_face = data_edge_firstparent_face[ t_e5 ];
      while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ] != t 
             &&
             data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ] != nullindex )
        current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ];
      assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ] != nullindex );
      assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ] == t );
      data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ] = t5; // FIXME
    }
    
    /* vertex t_v0: nothing needs to change */
    
    /* vertex t_v1: nothing needs to change */
    
    /* vertex t_v2: nothing needs to change */
    
    /* vertex t_v3: relink */ 
    if( data_vertex_firstparent_face[ t_v3 ] == t ) 
      data_vertex_firstparent_face[ t_v3 ] = t1; // FIXME
    else {
      int current_face = data_vertex_firstparent_face[ t_v2 ];
      while( data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ] != t 
             &&
             data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ] != nullindex )
        current_face = data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ];
      assert( data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ] != nullindex );
      assert( data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ] == t );
      data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ] = t1;
    }
    
    
    /* face f0: link */
    data_face_firstparent_tetrahedron[ f0 ] = t0; // FIXME
    
    /* face f1: link */
    data_face_firstparent_tetrahedron[ f1 ] = t0; // FIXME
    
    /* face f2: link */
    data_face_firstparent_tetrahedron[ f2 ] = t0; // FIXME
    
    /* face f3: link */
    data_face_firstparent_tetrahedron[ f3 ] = t0; // FIXME
    
    /* face f4: link */
    data_face_firstparent_tetrahedron[ f4 ] = t0; // FIXME
    
    /* face f5: link */
    data_face_firstparent_tetrahedron[ f5 ] = t0; // FIXME
    
    
    /* edge e0: link */
    data_edge_firstparent_tetrahedron[ e0 ] = t0; // FIXME
    data_edge_firstparent_face[ e0 ] = t0; // FIXME
    
    /* edge e1: link */
    data_edge_firstparent_tetrahedron[ e1 ] = t0; // FIXME
    data_edge_firstparent_face[ e1 ] = t0; // FIXME
    
    /* edge e2: link */
    data_edge_firstparent_tetrahedron[ e2 ] = t0; // FIXME
    data_edge_firstparent_face[ e2 ] = t1; // FIXME
    
    /* edge e3: link */
    data_edge_firstparent_tetrahedron[ e3 ] = t0; // FIXME
    data_edge_firstparent_face[ e3 ] = t1; // FIXME
    
    
    /* vertex vn: link */
    data_vertex_firstparent_tetrahedron[ counter_vertices ] = t0;
    data_vertex_firstparent_face       [ counter_vertices ] = f0;
    data_vertex_firstparent_edge       [ counter_vertices ] = e0;
    
    
    /* Counters */
    counter_tetrahedra = counter_tetrahedra + 3;
    counter_faces      = counter_faces      + 6;
    counter_edges      = counter_edges      + 4;
    counter_vertices   = counter_vertices   + 1;
    
    /* DONE */
    
    check();
}


void MeshSimplicial3D::midpoint_refinement_global()
{
    check();
    
    int N = counter_tetrahedra;
    
    for( int t = 0; t < N; t++ ) {
      
      midpoint_refinement( t );
      
    }
    
    check();
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




