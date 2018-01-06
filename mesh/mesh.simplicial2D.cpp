
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

#include "mesh.simplicial2D.hpp"




MeshSimplicial2D::MeshSimplicial2D( int outerdim )
:
    Mesh( 2, outerdim ),
    
    counter_triangles(0),
    counter_edges(0),
    counter_vertices(0),
    
    data_triangle_edges(0),
    data_edge_firstparent_triangle(0),
    data_triangle_nextparents_of_edges(0),
    
    data_triangle_vertices(0),
    data_vertex_firstparent_triangle(0),
    data_triangle_nextparents_of_vertices(0),
    
    data_edge_vertices(0),
    data_vertex_firstparent_edge(0),
    data_edge_nextparents_of_vertices(0)
{
    check();
}


// TODO: Update this function

MeshSimplicial2D::MeshSimplicial2D( 
    int outerdim,
    const Coordinates& coords,
    const std::vector<std::array<int,3>> triangle_vertices
)
:
    Mesh( 2, outerdim ),
    
    counter_triangles( triangle_vertices.size() ),
    counter_edges( 0 ),
    counter_vertices( 0 ),
    
    data_triangle_edges( triangle_vertices.size(), { nullindex, nullindex, nullindex } ),
    data_edge_firstparent_triangle( 0 ),
    data_triangle_nextparents_of_edges( triangle_vertices.size(), { nullindex, nullindex, nullindex } ),
    
    data_triangle_vertices( triangle_vertices ),
    data_vertex_firstparent_triangle( 0 ),
    data_triangle_nextparents_of_vertices( triangle_vertices.size(), { nullindex, nullindex, nullindex } ),
    
    data_edge_vertices( 0 ),
    data_vertex_firstparent_edge( 0 ),
    data_edge_nextparents_of_vertices( 0 )
{
    
    getcoordinates() = coords;
    
    /* 1. create all edges, allocate memory */
    
    data_edge_vertices.resize( counter_triangles * 3 );
    for( int t  = 0; t  < counter_triangles; t++  )
    for( int ei = 0; ei <                 3; ei++ )
    {
      data_edge_vertices[ 0 * counter_triangles + t ] = { data_triangle_vertices[t][0], data_triangle_vertices[t][1] };
      data_edge_vertices[ 1 * counter_triangles + t ] = { data_triangle_vertices[t][0], data_triangle_vertices[t][2] };
      data_edge_vertices[ 2 * counter_triangles + t ] = { data_triangle_vertices[t][1], data_triangle_vertices[t][2] };
    }
    
    std::sort( data_edge_vertices.begin(), data_edge_vertices.end() );
    auto it = std::unique( data_edge_vertices.begin(), data_edge_vertices.end() );
    data_edge_vertices.resize( it - data_edge_vertices.begin() );
    
    counter_edges = data_edge_vertices.size();
    
    data_edge_firstparent_triangle.resize   ( counter_edges, nullindex );
    data_edge_nextparents_of_vertices.resize( counter_edges, { nullindex, nullindex } );
    
    /* 2. Count vertices, allocate memory */
    
    counter_vertices = 0;
    for( const auto& duple : data_edge_vertices )
    for( const int& vertex : duple )
      counter_vertices = counter_vertices < vertex ? vertex : counter_vertices; 
    counter_vertices += 1;
    
    data_vertex_firstparent_triangle.resize( counter_vertices, nullindex );
    data_vertex_firstparent_edge.resize( counter_vertices, nullindex );
    
    for( auto t : data_triangle_vertices )
      std::cout << t[0] << space << t[1] << space << t[2] << std::endl;
    for( auto e : data_edge_vertices )
      std::cout << e[0] << space << e[1] << std::endl;
    std::cout << std::endl;
      
    /* 3. For each vertex, set the first parent triangle and the neighboring parent triangles */
    
    for( int t =  0; t  < counter_triangles; t++  )
    for( int vi = 0; vi <                 3; vi++ )
    {
      int v = data_triangle_vertices[t][vi];
      
      assert( 0 <= v && v < counter_vertices );
      
      if( data_vertex_firstparent_triangle[v] == nullindex ) {
        
        data_vertex_firstparent_triangle[v] = t;
        
      } else {
        
        int old_first_parent = data_vertex_firstparent_triangle[v];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_triangles );
        assert( data_triangle_nextparents_of_vertices[ t ][ vi ] == nullindex );
        
        data_vertex_firstparent_triangle[v] = t;
        data_triangle_nextparents_of_vertices[ t ][ vi ] = old_first_parent;
        
      }
      
      assert( data_vertex_firstparent_triangle[v] != nullindex );
      assert( 0 <= data_vertex_firstparent_triangle[v] && data_vertex_firstparent_triangle[v] < counter_triangles );  
    }
    
    /* 4. For each vertex, set the first parent edge and the neighboring parent edges */
    
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
    
    /* 5. For each edge, set the first parent triangle and the neighboring parent triangles */
    
    for( int t = 0; t < counter_triangles; t++ )
    for( int e = 0; e <     counter_edges; e++ )
    {
      int voe0 = data_edge_vertices[e][0];
      int voe1 = data_edge_vertices[e][1];
      
      int vot0 = data_triangle_vertices[t][0];
      int vot1 = data_triangle_vertices[t][1];
      int vot2 = data_triangle_vertices[t][2];
      
      int eot = nullindex;
      
      if( voe0 == vot0 && voe1 == vot1 ) eot = 0;
      if( voe0 == vot0 && voe1 == vot2 ) eot = 1;
      if( voe0 == vot1 && voe1 == vot2 ) eot = 2;
      
      if( eot == nullindex ) continue;
      
      data_triangle_edges[t][eot] = e;
      
      if( data_edge_firstparent_triangle[e] == nullindex ) {
        
        data_edge_firstparent_triangle[e] = t;
        
      } else {
        
        int old_first_parent = data_edge_firstparent_triangle[e];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_triangles );
        
        data_edge_firstparent_triangle[e] = t;
        
        assert( data_triangle_nextparents_of_edges[ t ][ eot ] == nullindex );
        data_triangle_nextparents_of_edges[ t ][ eot ] = old_first_parent;
        
      }
      
      assert( data_edge_firstparent_triangle[e] != nullindex );
      assert( 0 <= data_edge_firstparent_triangle[e] && data_edge_firstparent_triangle[e] < counter_triangles );  
    }
    
    
    
    
    check();
}


MeshSimplicial2D::MeshSimplicial2D( 
    int outerdim,
    const Coordinates& coords,
    const std::vector<std::array<int,3>> triangle_edges,
    const std::vector<int              > edge_firstparent_triangle,
    const std::vector<std::array<int,3>> triangle_nextparents_of_edges,
    const std::vector<std::array<int,3>> triangle_vertices,
    const std::vector<int              > vertex_firstparent_triangle,
    const std::vector<std::array<int,3>> triangle_nextparents_of_vertices,
    const std::vector<std::array<int,2>> edge_vertices,
    const std::vector<int              > vertex_firstparent_edge,
    const std::vector<std::array<int,2>> edge_nextparents_of_vertices    
)
:
    Mesh( 2, outerdim ),
    
    counter_triangles( triangle_vertices.size() ),
    counter_edges( edge_vertices.size() ),
    counter_vertices( vertex_firstparent_triangle.size() ),
    
    data_triangle_edges( triangle_edges ),
    data_edge_firstparent_triangle( edge_firstparent_triangle ),
    data_triangle_nextparents_of_edges( triangle_nextparents_of_edges ),
    
    data_triangle_vertices( triangle_vertices ),
    data_vertex_firstparent_triangle( vertex_firstparent_triangle ),
    data_triangle_nextparents_of_vertices( triangle_nextparents_of_vertices ),
    
    data_edge_vertices( edge_vertices ),
    data_vertex_firstparent_edge( vertex_firstparent_edge ),
    data_edge_nextparents_of_vertices( edge_nextparents_of_vertices )
{
    
    getcoordinates() = coords;
    
    check();
}


MeshSimplicial2D::~MeshSimplicial2D()
{
    
}


bool MeshSimplicial2D::operator== ( const MeshSimplicial2D& mesh ) const 
{
  return counter_triangles == mesh.counter_triangles
         &&
         counter_edges == mesh.counter_edges
         &&
         counter_vertices == mesh.counter_vertices
         &&
         data_triangle_edges == mesh.data_triangle_edges
         &&
         data_edge_firstparent_triangle == mesh.data_edge_firstparent_triangle
         &&
         data_triangle_nextparents_of_edges == mesh.data_triangle_nextparents_of_edges
         &&
         data_triangle_vertices == mesh.data_triangle_vertices 
         &&
         data_vertex_firstparent_triangle == mesh.data_vertex_firstparent_triangle
         &&
         data_triangle_nextparents_of_vertices == mesh.data_triangle_nextparents_of_vertices
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

bool MeshSimplicial2D::operator!= ( const MeshSimplicial2D& mesh ) const 
{
  return ! ( *this == mesh );
}

void MeshSimplicial2D::check() const
{
    
    /* 1. Check the array sizes */
    
    assert( counter_triangles == data_triangle_edges.size() );
    assert( counter_triangles == data_triangle_nextparents_of_edges.size() );
    assert( counter_triangles == data_triangle_vertices.size() );
    assert( counter_triangles == data_triangle_nextparents_of_vertices.size() );
    assert( counter_edges == data_edge_vertices.size() );
    assert( counter_edges == data_edge_nextparents_of_vertices.size() );
    assert( counter_edges == data_edge_firstparent_triangle.size() );
    assert( counter_vertices == data_vertex_firstparent_edge.size() );
    assert( counter_vertices == data_vertex_firstparent_triangle.size() );
    
    assert( count_vertices() == getcoordinates().getnumber() );
    
    
    /* 
     * each triangle: each edge is a valid index
     * each triangle: each edge is unique 
     * each triangle: the next parents are unique
     * each triangle: the next parents are actually parents 
     */
    
    for( int t = 0; t < counter_triangles; t++ )
    {
        
        assert( data_triangle_edges[t][0] != nullindex );
        assert( data_triangle_edges[t][1] != nullindex );
        assert( data_triangle_edges[t][2] != nullindex );
        
        assert( 0 <= data_triangle_edges[t][0] && data_triangle_edges[t][0] < counter_edges );
        assert( 0 <= data_triangle_edges[t][1] && data_triangle_edges[t][1] < counter_edges );
        assert( 0 <= data_triangle_edges[t][2] && data_triangle_edges[t][2] < counter_edges );
        
        assert( data_triangle_edges[t][0] != data_triangle_edges[t][1] );
        assert( data_triangle_edges[t][0] != data_triangle_edges[t][2] );
        assert( data_triangle_edges[t][1] != data_triangle_edges[t][2] );
        
        
        if( data_triangle_nextparents_of_edges[t][0] != nullindex || data_triangle_nextparents_of_edges[t][1] != nullindex )
          assert( data_triangle_nextparents_of_edges[t][0] != data_triangle_nextparents_of_edges[t][1] );
        
        if( data_triangle_nextparents_of_edges[t][0] != nullindex || data_triangle_nextparents_of_edges[t][2] != nullindex )
          assert( data_triangle_nextparents_of_edges[t][0] != data_triangle_nextparents_of_edges[t][2] );
        
        if( data_triangle_nextparents_of_edges[t][1] != nullindex || data_triangle_nextparents_of_edges[t][2] != nullindex )
          assert( data_triangle_nextparents_of_edges[t][1] != data_triangle_nextparents_of_edges[t][2] );
        
        for( int ei = 0; ei < 3; ei++ )
        {
            
            if( data_triangle_nextparents_of_edges[t][ei] != nullindex )
              assert( 0 <= data_triangle_nextparents_of_edges[t][ei] && data_triangle_nextparents_of_edges[t][ei] < counter_triangles );
        
            if( data_triangle_nextparents_of_edges[t][ei] != nullindex )
              assert( data_triangle_edges[ data_triangle_nextparents_of_edges[t][ei] ][0] == data_triangle_edges[t][ei] 
                      ||
                      data_triangle_edges[ data_triangle_nextparents_of_edges[t][ei] ][1] == data_triangle_edges[t][ei] 
                      ||
                      data_triangle_edges[ data_triangle_nextparents_of_edges[t][ei] ][2] == data_triangle_edges[t][ei] );
            
        }
        
        
    }
    
    
    /*
     * each triangle: each vertex is a valid index
     * each triangle: each vertex is unique 
     * each triangle: the next parents are unique
     * each triangle: the next parents are actually parents 
     */
    
    for( int t = 0; t < counter_triangles; t++ )
    {
        
        assert( data_triangle_vertices[t][0] != nullindex );
        assert( data_triangle_vertices[t][1] != nullindex );
        assert( data_triangle_vertices[t][2] != nullindex );
        
        assert( 0 <= data_triangle_vertices[t][0] && data_triangle_vertices[t][0] < counter_vertices );
        assert( 0 <= data_triangle_vertices[t][1] && data_triangle_vertices[t][1] < counter_vertices );
        assert( 0 <= data_triangle_vertices[t][2] && data_triangle_vertices[t][2] < counter_vertices );
        
        assert( data_triangle_vertices[t][0] != data_triangle_vertices[t][1] );
        assert( data_triangle_vertices[t][0] != data_triangle_vertices[t][2] );
        assert( data_triangle_vertices[t][1] != data_triangle_vertices[t][2] );
        
        
//         if( data_triangle_nextparents_of_vertices[t][0] != nullindex || data_triangle_nextparents_of_vertices[t][1] != nullindex )
//           assert( data_triangle_nextparents_of_vertices[t][0] != data_triangle_nextparents_of_vertices[t][1] );
//         
//         if( data_triangle_nextparents_of_vertices[t][0] != nullindex || data_triangle_nextparents_of_vertices[t][2] != nullindex )
//           assert( data_triangle_nextparents_of_vertices[t][0] != data_triangle_nextparents_of_vertices[t][2] );
//         
//         if( data_triangle_nextparents_of_vertices[t][1] != nullindex || data_triangle_nextparents_of_vertices[t][2] != nullindex )
//           assert( data_triangle_nextparents_of_vertices[t][1] != data_triangle_nextparents_of_vertices[t][2] );
        
        
        for( int vi = 0; vi < 3; vi++ )
        {
            
            if( data_triangle_nextparents_of_vertices[t][vi] != nullindex )
              assert( 0 <= data_triangle_nextparents_of_vertices[t][vi] && data_triangle_nextparents_of_vertices[t][vi] < counter_triangles );
        
            if( data_triangle_nextparents_of_vertices[t][vi] != nullindex )
              assert( data_triangle_vertices[ data_triangle_nextparents_of_vertices[t][vi] ][0] == data_triangle_vertices[t][vi] 
                      ||
                      data_triangle_vertices[ data_triangle_nextparents_of_vertices[t][vi] ][1] == data_triangle_vertices[t][vi] 
                      ||
                      data_triangle_vertices[ data_triangle_nextparents_of_vertices[t][vi] ][2] == data_triangle_vertices[t][vi] );
            
        }
        
        
    }
    
    
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
     * check that all triangles are unique, even up to permutation
     * check both the vertices and the edges listed for each triangle
     */
    
    for( int t1 = 0; t1 < counter_triangles; t1++ )
    for( int t2 = 0; t2 < counter_triangles; t2++ )
    {
        if( t1 == t2 ) continue;
        
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][0] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][1] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][2] );
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][0] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][2] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][1] );
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][1] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][0] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][2] );
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][1] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][2] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][0] );
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][2] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][0] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][1] );
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][2] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][1] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][0] );
        
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][0] || data_triangle_edges[t1][1] != data_triangle_edges[t2][1] || data_triangle_edges[t1][2] != data_triangle_edges[t2][2] );
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][0] || data_triangle_edges[t1][1] != data_triangle_edges[t2][2] || data_triangle_edges[t1][2] != data_triangle_edges[t2][1] );
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][1] || data_triangle_edges[t1][1] != data_triangle_edges[t2][0] || data_triangle_edges[t1][2] != data_triangle_edges[t2][2] );
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][1] || data_triangle_edges[t1][1] != data_triangle_edges[t2][2] || data_triangle_edges[t1][2] != data_triangle_edges[t2][0] );
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][2] || data_triangle_edges[t1][1] != data_triangle_edges[t2][0] || data_triangle_edges[t1][2] != data_triangle_edges[t2][1] );
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][2] || data_triangle_edges[t1][1] != data_triangle_edges[t2][1] || data_triangle_edges[t1][2] != data_triangle_edges[t2][0] );
        
        
    }
    
    
    
    
    
    
    /*
     * each triangle: each edge is listed correctly
     */
    
    for( int t = 0; t < counter_triangles; t++ )
    {
        
        assert( data_edge_vertices[ data_triangle_edges[t][0] ][0] == data_triangle_vertices[t][0] );
        assert( data_edge_vertices[ data_triangle_edges[t][0] ][1] == data_triangle_vertices[t][1] );
        assert( data_edge_vertices[ data_triangle_edges[t][1] ][0] == data_triangle_vertices[t][0] );
        assert( data_edge_vertices[ data_triangle_edges[t][1] ][1] == data_triangle_vertices[t][2] );
        assert( data_edge_vertices[ data_triangle_edges[t][2] ][0] == data_triangle_vertices[t][1] );
        assert( data_edge_vertices[ data_triangle_edges[t][2] ][1] == data_triangle_vertices[t][2] );
        
    }
    
    
    
    
    
    
    
    /* 
     * each first parent triangle of an edge: first parent is non-null and a valid parent
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
        int p = data_edge_firstparent_triangle[e];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_triangles );
        
        assert( data_triangle_edges[p][0] == e || data_triangle_edges[p][1] == e || data_triangle_edges[p][2] == e );
    }
    
    /* 
     * each first parent triangle of a vertex: first parent is non-null and a valid parent
     */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
        int p = data_vertex_firstparent_triangle[v];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_triangles );
        
        assert( data_triangle_vertices[p][0] == v || data_triangle_vertices[p][1] == v || data_triangle_vertices[p][2] == v );
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
    
    
    
    
    /*
     * check that each parent triangle of an edge is listed as a parent 
     */
    
    for( int t  = 0; t  < counter_triangles; t++ )
    for( int ei = 0; ei <                 3; ei++ )
    {
      
      int e = data_triangle_edges[t][ei];
      
      int p = data_edge_firstparent_triangle[e];
      
      assert( p != nullindex );
      
      while( p != t && p != nullindex )
        if( data_triangle_edges[p][0] == e )
          p = data_triangle_nextparents_of_edges[p][0];
        else if( data_triangle_edges[p][1] == e )
          p = data_triangle_nextparents_of_edges[p][1];
        else if( data_triangle_edges[p][2] == e )
          p = data_triangle_nextparents_of_edges[p][2];
        else
          assert(false);
        
      assert( p == t );
      
    }
    
    /* 
     * check that each parent triangle of a vertex is listed as a parent 
     */
    
    for( int t  = 0; t  < counter_triangles; t++ )
    for( int vi = 0; vi <                 3; vi++ )
    {
      
      int v = data_triangle_vertices[t][vi];
      
      int p = data_vertex_firstparent_triangle[v];
      
      assert( p != nullindex );
      
      while( p != t && p != nullindex )
        if( data_triangle_vertices[p][0] == v )
          p = data_triangle_nextparents_of_vertices[p][0];
        else if( data_triangle_vertices[p][1] == v )
          p = data_triangle_nextparents_of_vertices[p][1];
        else if( data_triangle_vertices[p][2] == v )
          p = data_triangle_nextparents_of_vertices[p][2];
        else
          assert(false);
        
      assert( p == t );
      
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






void MeshSimplicial2D::print( std::ostream& os ) const
{
    os << "Printe Triangulation of 2D Manifold!" << std::endl;
    
    os << counter_triangles << space << counter_edges << space << counter_vertices << nl;
    
    
    
    os << "Triangle edges" << std::endl;
    
    for( const auto& triple : data_triangle_edges )
      std::cout << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    os << "Edge first parent triangles" << std::endl;
    
    for( int fp : data_edge_firstparent_triangle )
      std::cout << fp << nl;
    
    os << "Triangle next parents of edges" << std::endl;
    
    for( const auto& triple : data_triangle_nextparents_of_edges )
      std::cout << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    
    
    
    os << "Triangle vertices" << std::endl;
    
    for( const auto& triple : data_triangle_vertices )
      std::cout << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    os << "Edge first parent triangles" << std::endl;
    
    for( int fp : data_vertex_firstparent_triangle )
      std::cout << fp << nl;
    
    os << "Triangle next parents of edges" << std::endl;
    
    for( const auto& triple : data_triangle_nextparents_of_vertices )
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






bool MeshSimplicial2D::dimension_counted( int dim ) const
{
    assert( 0 <= dim && dim <= 2 );
    return true;
}

int MeshSimplicial2D::count_simplices( int dim ) const
{
  if( dim == 0 )
    return count_vertices();
  else if( dim == 1 )
    return count_edges();
  else if( dim == 2 )
    return count_triangles();
  else
    assert(false);
}

bool MeshSimplicial2D::subsimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub <= sup && sup <= 2 );
    return true;
}

IndexMap MeshSimplicial2D::getsubsimplices( int sup, int sub, int cell ) const
{
  
  if( sup == 2 && sub == 2 ) {
    
    assert( 0 <= cell && cell < count_triangles() );
    return IndexMap( IndexRange(0,0), IndexRange(0,count_triangles()-1), { cell } );
    
  } else if( sup == 2 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_triangles() );
    auto temp = get_triangle_edges(cell);
    return IndexMap( IndexRange(0,2), IndexRange( 0, count_edges()-1 ), std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 2 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_triangles() );
    auto temp = get_triangle_vertices(cell);
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

bool MeshSimplicial2D::supersimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub <= sup && sup <= 2 );
    return true;
}

const std::vector<int> MeshSimplicial2D::getsupersimplices( int sup, int sub, int cell ) const
{
  
  if( sup == 2 && sub == 2 ) {
    
    assert( 0 <= cell && cell < count_triangles() );
    return { cell };
    
  } else if( sup == 2 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    auto temp = get_triangle_parents_of_edge( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 2 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    auto temp = get_triangle_parents_of_vertex( cell ); 
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

int MeshSimplicial2D::count_triangles() const
{
    return counter_triangles;
}

int MeshSimplicial2D::count_edges() const
{
    return counter_edges;
}

int MeshSimplicial2D::count_vertices() const
{
    return counter_vertices;
}




// TODO: Proofread the following methods 

/* subsimplex relation of triangles and edges */

bool MeshSimplicial2D::contains_triangle_edge( int t, int e ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= e && e < counter_edges );
    
    return ( data_triangle_edges[t][0] == e ) || ( data_triangle_edges[t][1] == e ) || ( data_triangle_edges[t][2] == e );
} 

int MeshSimplicial2D::indexof_triangle_edge( int t, int e ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= e && e < counter_edges );
    if     ( data_triangle_edges[t][0] == e ) return 0;
    else if( data_triangle_edges[t][1] == e ) return 1;
    else if( data_triangle_edges[t][2] == e ) return 2;
    else                                      assert(false);
} 

int MeshSimplicial2D::get_triangle_edge( int t, int ei ) const
{
    assert( 0 <= t  && t  < counter_triangles );
    assert( 0 <= ei && ei < 3 );
    return data_triangle_edges[t][ei];
} 

const std::array<int,3> MeshSimplicial2D::get_triangle_edges( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    return data_triangle_edges[t];
} 



/* subsimplex relation of triangle and vertices */

bool MeshSimplicial2D::contains_triangle_vertex( int t, int v ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= v && v < counter_vertices );
    
    return ( data_triangle_vertices[t][0] == v ) || ( data_triangle_vertices[t][1] == v ) || ( data_triangle_vertices[t][2] == v );
} 

int MeshSimplicial2D::indexof_triangle_vertex( int t, int v ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= v && v < counter_vertices );
    if     ( data_triangle_vertices[t][0] == v ) return 0;
    else if( data_triangle_vertices[t][1] == v ) return 1;
    else if( data_triangle_vertices[t][2] == v ) return 2;
    else                                         assert(false);
} 

int MeshSimplicial2D::get_triangle_vertex( int t, int vi ) const
{
    assert( 0 <= t  && t  < counter_triangles );
    assert( 0 <= vi && vi < 3 );
    return data_triangle_vertices[t][vi];
} 

const std::array<int,3> MeshSimplicial2D::get_triangle_vertices( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    return data_triangle_vertices[t];
} 




/* subsimplex relation of edges and vertices */

bool MeshSimplicial2D::contains_edge_vertex( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    
    return ( data_edge_vertices[e][0] == v ) || ( data_edge_vertices[e][1] == v );
} 

int MeshSimplicial2D::indexof_edge_vertex( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    if     ( data_edge_vertices[e][0] == v ) return 0;
    else if( data_edge_vertices[e][1] == v ) return 1;
    else                                     assert(false);
} 

int MeshSimplicial2D::get_edge_vertex( int e, int vi ) const
{
    assert( 0 <= e  && e  < counter_edges );
    assert( 0 <= vi && vi < 2 );
    return data_edge_vertices[e][vi];
} 

const std::array<int,2> MeshSimplicial2D::get_edge_vertices( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    return data_edge_vertices[e];
} 




// TODO: Proofread the following methods 

/* triangle parents of a edge */

int MeshSimplicial2D::count_edge_triangle_parents( int e ) const
{
  return get_triangle_parents_of_edge( e ).size();
}

int MeshSimplicial2D::get_edge_firstparent_triangle( int e ) const
{
  assert( 0 <= e && e < counter_edges );
  return data_edge_firstparent_triangle[ e ];
}

int MeshSimplicial2D::get_edge_nextparent_triangle( int e, int t ) const
{
  assert( 0 <= e && e < counter_edges );
  assert( 0 <= t && t < counter_triangles );
  
  if( data_triangle_edges[t][0] == e )
    return data_triangle_nextparents_of_edges[t][0];
  else if( data_triangle_edges[t][1] == e )
    return data_triangle_nextparents_of_edges[t][1];
  else if( data_triangle_edges[t][2] == e )
    return data_triangle_nextparents_of_edges[t][2];
  else
    assert(false);
}

int MeshSimplicial2D::get_triangle_nextparent_of_edge( int t, int ei ) const
{
  assert( 0 <= t  && t  < counter_triangles );
  assert( 0 <= ei && ei < 3 );
  return data_triangle_nextparents_of_edges[t][ei];
}

bool MeshSimplicial2D::is_triangle_edge_parent( int t, int e ) const
{
  assert( 0 <= e && e < counter_edges );
  assert( 0 <= t && t < counter_triangles );
  return data_triangle_edges[t][0] == e || data_triangle_edges[t][1] == e || data_triangle_edges[t][2] == e;
}

int MeshSimplicial2D::indexof_triangle_edge_parent( int t, int e ) const
{
  assert( 0 <= e && e < count_edges() );
  std::vector<int> triangles = get_triangle_parents_of_edge( e );
  
  auto iter = std::find( triangles.begin(), triangles.end(), t ); 
  assert( iter != triangles.end() );
  
  return iter - triangles.begin();
}

std::vector<int> MeshSimplicial2D::get_triangle_parents_of_edge( int e ) const
{
  assert( 0 <= e && e < count_edges() );
  
  std::vector<int> ret;
  
  for( int t = 0; t < count_triangles(); t++ ) 
    for( int te : get_triangle_edges(t) )
      if( e == te )
        ret.push_back( t );
  
  return ret;
}


/* triangle parents of a vertex */

int MeshSimplicial2D::count_vertex_triangle_parents( int v ) const
{
  return get_triangle_parents_of_vertex( v ).size();
}

int MeshSimplicial2D::get_vertex_firstparent_triangle( int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  return data_vertex_firstparent_triangle[ v ];
}

int MeshSimplicial2D::get_vertex_nextparent_triangle( int v, int t ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= t && t < counter_triangles );
  
  if( data_triangle_vertices[t][0] == v )
    return data_triangle_nextparents_of_vertices[t][0];
  else if( data_triangle_vertices[t][1] == v )
    return data_triangle_nextparents_of_vertices[t][1];
  else if( data_triangle_vertices[t][2] == v )
    return data_triangle_nextparents_of_vertices[t][2];
  else
    assert(false);
}

int MeshSimplicial2D::get_triangle_nextparent_of_vertex( int t, int vi ) const
{
  assert( 0 <= t  && t  < counter_triangles );
  assert( 0 <= vi && vi < 3 );
  return data_triangle_nextparents_of_vertices[t][vi];
}

bool MeshSimplicial2D::is_triangle_vertex_parent( int t, int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= t && t < counter_triangles );
  return data_triangle_vertices[t][0] == v || data_triangle_vertices[t][1] == v || data_triangle_vertices[t][2] == v;
}

int MeshSimplicial2D::indexof_triangle_vertex_parent( int t, int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  std::vector<int> triangles = get_triangle_parents_of_vertex( v );
  
  auto iter = std::find( triangles.begin(), triangles.end(), t ); 
  assert( iter != triangles.end() );
  
  return iter - triangles.begin();
}

std::vector<int> MeshSimplicial2D::get_triangle_parents_of_vertex( int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  
  std::vector<int> ret;
  
  for( int t = 0; t < count_triangles(); t++ ) 
    for( int tv : get_triangle_vertices(t) )
      if( v == tv )
        ret.push_back( t );
  
  return ret;
}




/* edge parents of a vertex */

int MeshSimplicial2D::count_vertex_edge_parents( int v ) const
{
  return get_edge_parents_of_vertex( v ).size();
}

int MeshSimplicial2D::get_vertex_firstparent_edge( int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  return data_vertex_firstparent_edge[ v ];
}

int MeshSimplicial2D::get_vertex_nextparent_edge( int v, int e ) const
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

int MeshSimplicial2D::get_edge_nextparent_of_vertex( int e, int vi ) const
{
  assert( 0 <= e  && e  < counter_edges );
  assert( 0 <= vi && vi < 2 );
  return data_edge_nextparents_of_vertices[e][vi];
}

bool MeshSimplicial2D::is_edge_vertex_parent( int e, int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= e && e < counter_edges    );
  return data_edge_vertices[e][0] == v || data_edge_vertices[e][1] == v;
}

int MeshSimplicial2D::indexof_edge_vertex_parent( int e, int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  std::vector<int> edges = get_edge_parents_of_vertex( v );
  
  auto iter = std::find( edges.begin(), edges.end(), e ); 
  assert( iter != edges.end() );
  
  return iter - edges.begin();
}

std::vector<int> MeshSimplicial2D::get_edge_parents_of_vertex( int v ) const
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
 * * * * * NEWEST VERTEX BISECTION
 */
    
void MeshSimplicial2D::bisect_edge( int e )
{
    assert( 0 <= e && e < counter_edges );
    check();
    
    /* 
     * DESCRIPTION
     * 
     * The Bisection works as follows:
     * 
     * triangle -> edge 
     * 
     *          T2E: 
     *          
     *            just fill in the data. 
     *          
     *          next parents + first parent:
     *            
     *            for outer edges:
     *            start with the first parent pointer and run over all parents,
     *            replacing the old triangle with the new triangle if need be.
     *            
     *            for inner edges: 
     *            fill in the data
     *            
     *            Finally, create parent data for the two new edges 
     *          
     * 
     * triangle -> vertex 
     *          
     *          T2V: 
     *          
     *            just fill in the data.
     *          
     *          next parents + first parent:
     *          
     *            for front,back,opposing vertex:
     *            start with the first parent pointer and run over all parents,
     *            replacing the old triangle with the new triangle if need be.
     *            
     *            Finally, create new data for new vertex.
     *          
     *  
     * edge -> vertex 
     * 
     *          E2V: 
     *          
     *            just fill in the data.
     *            
     *          next parents + first parent:
     *          
     *            for back vertex:
     *                  [nothing changes]
     *            for front vertex:
     *                  run over all parents and replace the old edge.
     *            for opposing vertex:
     *                  add the additional edge as first parent
     * 
     * 
     * */
    
    
    
    
    /*
     * DATA COLLECTION 
     */
    
    int e_back_vertex  = data_edge_vertices[ e ][ 0 ];
    int e_front_vertex = data_edge_vertices[ e ][ 1 ];
    
    std::vector<int> old_triangles = get_triangle_parents_of_edge( e );
    
    std::vector<int> localindex_of_refinementedge( old_triangles.size() );
    
    
    FloatVector midcoordinate = get_edge_midpoint( e );
    
    
    
    /*
     * ALLOCATE MEMORY FOR THE DATA  
     */
    
    data_triangle_nextparents_of_edges.resize( counter_triangles + old_triangles.size()    , { nullindex, nullindex, nullindex } );
    data_triangle_edges.resize               ( counter_triangles + old_triangles.size()    , { nullindex, nullindex, nullindex } );
    data_edge_firstparent_triangle.resize    ( counter_edges     + old_triangles.size() + 1,                           nullindex );
    
    data_triangle_nextparents_of_vertices.resize( counter_triangles + old_triangles.size(), { nullindex, nullindex, nullindex } );
    data_triangle_vertices.resize               ( counter_triangles + old_triangles.size(), { nullindex, nullindex, nullindex } );
    data_vertex_firstparent_triangle.resize     ( counter_vertices  + 1,                                              nullindex );
    
    data_edge_nextparents_of_vertices.resize( counter_edges    + old_triangles.size() + 1, { nullindex, nullindex } );
    data_edge_vertices.resize               ( counter_edges    + old_triangles.size() + 1, { nullindex, nullindex } );
    data_vertex_firstparent_edge.resize     ( counter_vertices + 1,                                       nullindex );
    
    
    
    
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
    
    // no parent triangles yet for the new vertex; will be filled in below */
    data_vertex_firstparent_triangle[ counter_vertices ] = nullindex;
    
    // no parent triangles yet for the front edge; will be filled in below */
    data_edge_firstparent_triangle[ counter_edges ] = nullindex;
    
    
    
    std::cout << "BIG LOOP" << std::endl;
    
    for( int ot = 0; ot < old_triangles.size(); ot++ ) {
      
      int t_old = old_triangles[ ot ];
      int t_new = counter_triangles + ot;
      
      int t_e0 = data_triangle_edges[ t_old ][ 0 ];
      int t_e1 = data_triangle_edges[ t_old ][ 1 ];
      int t_e2 = data_triangle_edges[ t_old ][ 2 ];
      
      int t_v0 = data_triangle_vertices[ t_old ][ 0 ];
      int t_v1 = data_triangle_vertices[ t_old ][ 1 ];
      int t_v2 = data_triangle_vertices[ t_old ][ 2 ];
      
      int t_e_n0 = data_triangle_nextparents_of_edges[ t_old ][ 0 ];
      int t_e_n1 = data_triangle_nextparents_of_edges[ t_old ][ 1 ];
      int t_e_n2 = data_triangle_nextparents_of_edges[ t_old ][ 2 ];
      
      int t_v_n0 = data_triangle_nextparents_of_vertices[ t_old ][ 0 ];
      int t_v_n1 = data_triangle_nextparents_of_vertices[ t_old ][ 1 ];
      int t_v_n2 = data_triangle_nextparents_of_vertices[ t_old ][ 2 ];
      
      localindex_of_refinementedge[ ot ] = ( t_e0 == e ) ? 0 : ( ( t_e1 == e ) ? 1 : 2 );
      assert( data_triangle_edges[ t_old ][ localindex_of_refinementedge[ ot ] ] == e );
      
      
      
      if( localindex_of_refinementedge[ ot ] == 0 ) { // 0 1 
        
        assert( t_v0 == e_back_vertex && t_v1 == e_front_vertex && e == t_e0 );
        
        /* triangle vertices */
        data_triangle_vertices[ t_old ][0] = t_v0;
        data_triangle_vertices[ t_old ][1] = counter_vertices;
        data_triangle_vertices[ t_old ][2] = t_v2;
        
        data_triangle_vertices[ t_new ][0] = counter_vertices;
        data_triangle_vertices[ t_new ][1] = t_v1;
        data_triangle_vertices[ t_new ][2] = t_v2;
        
        /* triangle edges */
        data_triangle_edges[ t_old ][0] = e;
        data_triangle_edges[ t_old ][1] = t_e1;
        data_triangle_edges[ t_old ][2] = counter_edges + 1 + ot;
        
        data_triangle_edges[ t_new ][0] = counter_edges;
        data_triangle_edges[ t_new ][1] = counter_edges + 1 + ot;
        data_triangle_edges[ t_new ][2] = t_e2;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + ot ][0] = counter_vertices;
        data_edge_vertices[ counter_edges + 1 + ot ][1] = t_v2;
        
        
        /* triangle vertex neighbors */
        data_triangle_nextparents_of_vertices[ t_old ][0] = t_v_n0;
        data_triangle_nextparents_of_vertices[ t_old ][1] = counter_triangles + ot;
        data_triangle_nextparents_of_vertices[ t_old ][2] = counter_triangles + ot;
        
        data_triangle_nextparents_of_vertices[ t_new ][0] = data_vertex_firstparent_triangle[ counter_vertices ];
        data_triangle_nextparents_of_vertices[ t_new ][1] = t_v_n1;
        data_triangle_nextparents_of_vertices[ t_new ][2] = t_v_n2;
        
        data_vertex_firstparent_triangle[ counter_vertices ] = t_old;
        
        /* triangle edge neighbors */
        data_triangle_nextparents_of_edges[ t_old ][0] = t_e_n0;
        data_triangle_nextparents_of_edges[ t_old ][1] = t_e_n1;
        data_triangle_nextparents_of_edges[ t_old ][2] = counter_triangles + ot;
        
        data_triangle_nextparents_of_edges[ t_new ][0] = data_edge_firstparent_triangle[ counter_edges ];
        data_triangle_nextparents_of_edges[ t_new ][1] = nullindex;
        data_triangle_nextparents_of_edges[ t_new ][2] = t_e_n2;
        
        data_edge_firstparent_triangle[ counter_edges ] = t_new;
        
        data_edge_firstparent_triangle[ counter_edges + 1 + ot ] = t_old;
        
        /* run over the front outer edge parent triangle list and replace 't_old' by 't_new' */
        if( data_edge_firstparent_triangle[ t_e2 ] == t_old ) 
          data_edge_firstparent_triangle[ t_e2 ] = t_new;
        else {
          int current_triangle = data_edge_firstparent_triangle[ t_e2 ];
          while( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != t_old 
                 &&
                 data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != nullindex )
            current_triangle = data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ];
          assert( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != nullindex );
          data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] = t_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][0] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ t_v2 ];
        data_vertex_firstparent_edge[ t_v2 ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][1] = opposing_vertex_firstparent_edge;
        
        
        
      } else if( localindex_of_refinementedge[ ot ] == 1 ) { // 0 2 
        
        assert( t_v0 == e_back_vertex && t_v2 == e_front_vertex && e == t_e1 );
        
        /* triangle vertices */
        data_triangle_vertices[ t_old ][0] = t_v0;
        data_triangle_vertices[ t_old ][1] = t_v1;
        data_triangle_vertices[ t_old ][2] = counter_vertices;
        
        data_triangle_vertices[ t_new ][0] = t_v1;
        data_triangle_vertices[ t_new ][1] = counter_vertices;
        data_triangle_vertices[ t_new ][2] = t_v2;
        
        /* triangle edges */
        data_triangle_edges[ t_old ][0] = t_e0;
        data_triangle_edges[ t_old ][1] = e;
        data_triangle_edges[ t_old ][2] = counter_edges + 1 + ot;
        
        data_triangle_edges[ t_new ][0] = counter_edges + 1 + ot;
        data_triangle_edges[ t_new ][1] = t_e2;
        data_triangle_edges[ t_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + ot ][0] = t_v1;
        data_edge_vertices[ counter_edges + 1 + ot ][1] = counter_vertices;
        
        
        /* triangle vertex neighbors */
        data_triangle_nextparents_of_vertices[ t_old ][0] = t_v_n0;
        data_triangle_nextparents_of_vertices[ t_old ][1] = counter_triangles + ot;
        data_triangle_nextparents_of_vertices[ t_old ][2] = counter_triangles + ot;
        
        data_triangle_nextparents_of_vertices[ t_new ][0] = t_v_n1;
        data_triangle_nextparents_of_vertices[ t_new ][1] = data_vertex_firstparent_triangle[ counter_vertices ];
        data_triangle_nextparents_of_vertices[ t_new ][2] = t_v_n2;
        
        data_vertex_firstparent_triangle[ counter_vertices ] = t_old;
        
        /* triangle edge neighbors */
        data_triangle_nextparents_of_edges[ t_old ][0] = t_e_n0;
        data_triangle_nextparents_of_edges[ t_old ][1] = t_e_n1;
        data_triangle_nextparents_of_edges[ t_old ][2] = counter_triangles + ot;
        
        data_triangle_nextparents_of_edges[ t_new ][0] = nullindex;
        data_triangle_nextparents_of_edges[ t_new ][1] = t_e_n2;
        data_triangle_nextparents_of_edges[ t_new ][2] = data_edge_firstparent_triangle[ counter_edges ];
        
        data_edge_firstparent_triangle[ counter_edges ] = t_new;
        
        data_edge_firstparent_triangle[ counter_edges + 1 + ot ] = t_old;
        
        /* run over the front outer edge parent triangle list and replace 't_old' by 't_new' */
        if( data_edge_firstparent_triangle[ t_e2 ] == t_old ) 
          data_edge_firstparent_triangle[ t_e2 ] = t_new;
        else {
          int current_triangle = data_edge_firstparent_triangle[ t_e2 ];
          while( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != t_old 
                 &&
                 data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != nullindex )
            current_triangle = data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ];
          assert( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != nullindex );
          data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] = t_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][1] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ t_v1 ];
        data_vertex_firstparent_edge[ t_v1 ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][0] = opposing_vertex_firstparent_edge;
        
        
      } else if( localindex_of_refinementedge[ ot ] == 2 ) { // 1 2 
        
        assert( t_v1 == e_back_vertex && t_v2 == e_front_vertex && e == t_e2 );
        
        /* triangle vertices */
        data_triangle_vertices[ t_old ][0] = t_v0;
        data_triangle_vertices[ t_old ][1] = t_v1;
        data_triangle_vertices[ t_old ][2] = counter_vertices;
        
        data_triangle_vertices[ t_new ][0] = t_v0;
        data_triangle_vertices[ t_new ][1] = counter_vertices;
        data_triangle_vertices[ t_new ][2] = t_v2;
        
        /* triangle edges */
        data_triangle_edges[ t_old ][0] = t_e0;
        data_triangle_edges[ t_old ][1] = counter_edges + 1 + ot;
        data_triangle_edges[ t_old ][2] = e;
        
        data_triangle_edges[ t_new ][0] = counter_edges + 1 + ot;
        data_triangle_edges[ t_new ][1] = t_e1;
        data_triangle_edges[ t_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + ot ][0] = t_v0;
        data_edge_vertices[ counter_edges + 1 + ot ][1] = counter_vertices;
        
        
        /* triangle vertex neighbors */
        data_triangle_nextparents_of_vertices[ t_old ][0] = counter_triangles + ot;
        data_triangle_nextparents_of_vertices[ t_old ][1] = t_v_n1;
        data_triangle_nextparents_of_vertices[ t_old ][2] = counter_triangles + ot;
        
        data_triangle_nextparents_of_vertices[ t_new ][0] = t_v_n0;
        data_triangle_nextparents_of_vertices[ t_new ][1] = data_vertex_firstparent_triangle[ counter_vertices ];
        data_triangle_nextparents_of_vertices[ t_new ][2] = t_v_n2;
        
        data_vertex_firstparent_triangle[ counter_vertices ] = t_old;
        
        /* triangle edge neighbors */
        data_triangle_nextparents_of_edges[ t_old ][0] = t_e_n0;
        data_triangle_nextparents_of_edges[ t_old ][1] = counter_triangles + ot;
        data_triangle_nextparents_of_edges[ t_old ][2] = t_e_n2;
        
        data_triangle_nextparents_of_edges[ t_new ][0] = nullindex;
        data_triangle_nextparents_of_edges[ t_new ][1] = t_e_n1;
        data_triangle_nextparents_of_edges[ t_new ][2] = data_edge_firstparent_triangle[ counter_edges ];
        
        data_edge_firstparent_triangle[ counter_edges ] = t_new;
        
        data_edge_firstparent_triangle[ counter_edges + 1 + ot ] = t_old;
        
        /* run over the front outer edge parent triangle list and replace 't_old' by 't_new' */
        if( data_edge_firstparent_triangle[ t_e1 ] == t_old ) 
          data_edge_firstparent_triangle[ t_e1 ] = t_new;
        else {
          int current_triangle = data_edge_firstparent_triangle[ t_e1 ];
          while( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] != t_old 
                 &&
                 data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] != nullindex )
            current_triangle = data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ];
          assert( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] != nullindex );
          data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] = t_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][1] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ t_v0 ];
        data_vertex_firstparent_edge[ t_v0 ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][0] = opposing_vertex_firstparent_edge;
        
      } else {
        
        assert(false);
        
      } 
      
    }
    
    
    std::cout << "BIG LOOP END" << std::endl;
    
    
    /* Run over the front vertex' parent triangles and conduct manipulations */

    int* pointer_to_index = &data_vertex_firstparent_triangle[ e_front_vertex ];
    
    while( *pointer_to_index != nullindex ) { 
      
      std::vector<int>::iterator it = std::find( old_triangles.begin(), old_triangles.end(), *pointer_to_index );
      
      if( it != old_triangles.end() ) {
        
        assert( *pointer_to_index == *it );
        assert( *it == old_triangles[ it - old_triangles.begin() ] );
        
        *pointer_to_index = counter_triangles + ( it - old_triangles.begin() );
        
        std::cout << "manipulate" << std::endl;
        
      } else std::cout << "keep" << std::endl;
      
      int localindex_of_front_vertex = nullindex;
      if( data_triangle_vertices[ *pointer_to_index ][ 0 ] == e_front_vertex ) localindex_of_front_vertex = 0;
      if( data_triangle_vertices[ *pointer_to_index ][ 1 ] == e_front_vertex ) localindex_of_front_vertex = 1;
      if( data_triangle_vertices[ *pointer_to_index ][ 2 ] == e_front_vertex ) localindex_of_front_vertex = 2;
      assert( localindex_of_front_vertex != nullindex );
      
      pointer_to_index = &( data_triangle_nextparents_of_vertices[ *pointer_to_index ][ localindex_of_front_vertex ] );
      
    }
    
    getcoordinates().append( midcoordinate );
    
    
    std::cout << "FINISHED" << std::endl;
    
    /*
     *  UPDATE COUNTERS 
     */
    
    counter_triangles += old_triangles.size();
    counter_edges     += 1 + old_triangles.size();
    counter_vertices  += 1;
    
    
    
    
    
    /* Done */
    
    check();
    
}

/*
 * * * * * NEWEST VERTEX BISECTION
 */
    







void MeshSimplicial2D::uniformrefinement()
{
    check();
    
    /* resize the arrays */
    
    data_triangle_edges.resize               ( 4 * counter_triangles                    , { nullindex, nullindex, nullindex } );
    data_edge_firstparent_triangle.resize    ( 2 * counter_edges + 3 * counter_triangles, nullindex                           );
    data_triangle_nextparents_of_edges.resize( 4 * counter_triangles                    , { nullindex, nullindex, nullindex } );
    
    data_triangle_vertices.resize               ( 4 * counter_triangles           , { nullindex, nullindex, nullindex } );
    data_vertex_firstparent_triangle.resize     ( counter_vertices + counter_edges, nullindex                           );
    data_triangle_nextparents_of_vertices.resize( 4 * counter_triangles           , { nullindex, nullindex, nullindex } );
    
    data_edge_vertices.resize               ( 2 * counter_edges + 3 * counter_triangles, { nullindex, nullindex } );
    data_vertex_firstparent_edge.resize     ( counter_vertices + counter_edges         , nullindex                );
    data_edge_nextparents_of_vertices.resize( 2 * counter_edges + 3 * counter_triangles, { nullindex, nullindex } );
    
    getcoordinates().addcoordinates( counter_edges );
    
    
    
    /* create the new coordinates and fill them up */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      getcoordinates().loadvector( counter_vertices + e, get_edge_midpoint( e ) );
    }
    
    
    
    
    /**** Edge -- Vertex Relation ****/
    
    /* for each old vertex, set the new parent edge */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      int p = data_vertex_firstparent_edge[v];
      
      assert( p != nullindex && 0 <= p && p < counter_edges );
      
      int vi = data_edge_vertices[p][0] == v ? 0 : 1;
      
      assert( data_edge_vertices[p][0] == v || data_edge_vertices[p][1] == v );
      assert( data_edge_vertices[p][vi] == v );
      
      data_vertex_firstparent_edge[v] = vi * counter_edges + p;
    }
    
    
    /* for each old edge, relocate the data of the old vertices' old parent edges */
    
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
    
    
    /* for each new vertex, 
     * set the first and second parent edge from the old edge 
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      data_vertex_firstparent_edge[counter_vertices + e] = e;
      
      data_edge_nextparents_of_vertices[ 0 * counter_edges + e ][1] = counter_edges + e;
      data_edge_nextparents_of_vertices[ 1 * counter_edges + e ][0] = nullindex;
    }
    
    
    /* for each old edge, run over the adjacent triangles 
     * and add the corresponding new edges to the list of 
     * parent edges of new vertex.
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int t = data_edge_firstparent_triangle[e];
      
      while( t != nullindex ) {
        
        int ei   = nullindex;
        int e_1  = nullindex; 
        int e_2  = nullindex;
        int vi_1 = nullindex;
        int vi_2 = nullindex;
        
        // [ 01 02 12 ] -> [ 01 02 ] [ 01 12 ] [ 02 12 ]
        
        if(        data_triangle_edges[t][0] == e ) {
          ei = 0; e_1 = 0; e_2 = 1; vi_1 = 0; vi_2 = 0; 
        } else if( data_triangle_edges[t][1] == e ) {
          ei = 1; e_1 = 0; e_2 = 2; vi_1 = 1; vi_2 = 0; 
        } else if( data_triangle_edges[t][2] == e ) {
          ei = 2; e_1 = 1; e_2 = 2; vi_1 = 1; vi_2 = 1; 
        } else
          assert(false);
        
        assert( ei  != nullindex && e_1 != nullindex && e_2 != nullindex );
        
        int old_first_parent = data_vertex_firstparent_edge[ counter_vertices + e ];
        
        data_vertex_firstparent_edge[ counter_vertices + e ]
          = 2 * counter_edges + e_1 * counter_triangles + t;
        
        data_edge_nextparents_of_vertices[ 2 * counter_edges + e_1 * counter_triangles + t ][ vi_1 ]
          = 2 * counter_edges + e_2 * counter_triangles + t;
        
        data_edge_nextparents_of_vertices[ 2 * counter_edges + e_2 * counter_triangles + t ][ vi_2 ]
          = old_first_parent;
        
        t = data_triangle_nextparents_of_edges[ t ][ ei ];
        
      }
      
    }
    
    /* for each edge created from an old edge, set the vertices */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int vertex_back  = data_edge_vertices[e][0];
      int vertex_front = data_edge_vertices[e][1];
      
      data_edge_vertices[e                ][0] = vertex_back;
      data_edge_vertices[e                ][1] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][0] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][1] = vertex_front;
    }
    
    /* for each triangle, set the vertices of the new edge */
    
    for( int t = 0; t < counter_triangles; t++ )
    {
      data_edge_vertices[ 2 * counter_edges + 0 * counter_triangles + t ][0] = counter_vertices + data_triangle_edges[t][0];
      data_edge_vertices[ 2 * counter_edges + 0 * counter_triangles + t ][1] = counter_vertices + data_triangle_edges[t][1];
      data_edge_vertices[ 2 * counter_edges + 1 * counter_triangles + t ][0] = counter_vertices + data_triangle_edges[t][0];
      data_edge_vertices[ 2 * counter_edges + 1 * counter_triangles + t ][1] = counter_vertices + data_triangle_edges[t][2];
      data_edge_vertices[ 2 * counter_edges + 2 * counter_triangles + t ][0] = counter_vertices + data_triangle_edges[t][1];
      data_edge_vertices[ 2 * counter_edges + 2 * counter_triangles + t ][1] = counter_vertices + data_triangle_edges[t][2];
    }
    
    
    /**** Triangle -- Vertex Relation ****/
    
    /* for each old vertex, set the new parent triangle */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      int p = data_vertex_firstparent_triangle[v];
      
      assert( p != nullindex && 0 <= p && p < counter_triangles );
      
      int vi = data_triangle_vertices[p][0] == v ? 0 : data_triangle_vertices[p][1] == v ? 1 : 2;
      
      assert( data_triangle_vertices[p][0] == v || data_triangle_vertices[p][1] == v  || data_triangle_vertices[p][2] == v );
      assert( data_triangle_vertices[p][vi] == v );
      
      data_vertex_firstparent_triangle[v] = vi * counter_triangles + p;
    }
    
    
    /* for each old triangle, relocate the data of the old vertices' parent triangle */
    
    for( int t  = 0; t  < counter_triangles;  t++ )
    for( int vi = 0; vi <                 3; vi++ )
    {
      int q = data_triangle_nextparents_of_vertices[t][vi];
      
      int v = data_triangle_vertices[t][vi];
      
      if( q == nullindex ) {
        
        data_triangle_nextparents_of_vertices[ vi * counter_triangles + t ][vi] = nullindex;
        
      } else {
        
        int vinp = data_triangle_vertices[q][0] == v ? 0 : data_triangle_vertices[q][1] == v ? 1 : 2;
        
        assert( data_triangle_vertices[q][0] == v || data_triangle_vertices[q][1] == v || data_triangle_vertices[q][2] == v );
        assert( data_triangle_vertices[q][vinp] == v );
        
        data_triangle_nextparents_of_vertices[ vi * counter_triangles + t ][vi] = vinp * counter_triangles + q;
      
      } 
      
    }
    
    /* for each old edge, run over the adjacent triangles 
     * and add the corresponding new triangles to the list of 
     * parent triangles of new vertex.
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int t = data_edge_firstparent_triangle[e];
      
      while( t != nullindex ) {
        
        int ei   = nullindex;
        int t_1  = nullindex; 
        int t_2  = nullindex;
        int t_3  = nullindex;
        int vi_1 = nullindex;
        int vi_2 = nullindex;
        int vi_3 = nullindex;
        
        // [ 01 02 12 ] -> [ 01 02 ] [ 01 12 ] [ 02 12 ]
        // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
        if(        data_triangle_edges[t][0] == e ) {
          ei = 0; 
          t_1 = 0; t_2 = 3; t_3 = 1; vi_1 = 1; vi_2 = 0; vi_3 = 0;  
        } else if( data_triangle_edges[t][1] == e ) {
          ei = 1; 
          t_1 = 0; t_2 = 3; t_3 = 2; vi_1 = 2; vi_2 = 1; vi_3 = 0;  
        } else if( data_triangle_edges[t][2] == e ) {
          ei = 2; 
          t_1 = 1; t_2 = 3; t_3 = 2; vi_1 = 2; vi_2 = 2; vi_3 = 1; 
        } else
          assert(false);
        
        int old_first_parent = data_vertex_firstparent_triangle[ counter_vertices + e ];
        
        data_vertex_firstparent_triangle[ counter_vertices + e ]
          = t_1 * counter_triangles + t;
        
        if( t_1 != 0 ) assert( data_triangle_nextparents_of_vertices[ t_1 * counter_triangles + t ][ vi_1 ] == nullindex );
        data_triangle_nextparents_of_vertices[ t_1 * counter_triangles + t ][ vi_1 ]
          = t_2 * counter_triangles + t;
        
        if( t_2 != 0 ) assert( data_triangle_nextparents_of_vertices[ t_2 * counter_triangles + t ][ vi_2 ] == nullindex );
        data_triangle_nextparents_of_vertices[ t_2 * counter_triangles + t ][ vi_2 ]
          = t_3 * counter_triangles + t;
        
        if( t_3 != 0 ) assert( data_triangle_nextparents_of_vertices[ t_3 * counter_triangles + t ][ vi_3 ] == nullindex );
        data_triangle_nextparents_of_vertices[ t_3 * counter_triangles + t ][ vi_3 ]
          = old_first_parent;
        
        t = data_triangle_nextparents_of_edges[ t ][ ei ];
        
      }
      
    }
    
    
    
//     /* for each triangle, create the new vertices */
//     
//     for( int t = 0; t < counter_triangles; t++ )
//     {
//       int v00 = data_triangle_vertices[t][0];
//       int v11 = data_triangle_vertices[t][1];
//       int v22 = data_triangle_vertices[t][2];
//       
//       int v01 = counter_vertices + data_triangle_edges[t][0];
//       int v02 = counter_vertices + data_triangle_edges[t][1];
//       int v12 = counter_vertices + data_triangle_edges[t][2];
//       
//       // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
//         
//       data_triangle_vertices[ 0 * counter_triangles + t ][0] = v00;
//       data_triangle_vertices[ 0 * counter_triangles + t ][1] = v01;
//       data_triangle_vertices[ 0 * counter_triangles + t ][2] = v02;
//       
//       data_triangle_vertices[ 1 * counter_triangles + t ][0] = v01;
//       data_triangle_vertices[ 1 * counter_triangles + t ][1] = v11;
//       data_triangle_vertices[ 1 * counter_triangles + t ][2] = v12;
//       
//       data_triangle_vertices[ 2 * counter_triangles + t ][0] = v02;
//       data_triangle_vertices[ 2 * counter_triangles + t ][1] = v12;
//       data_triangle_vertices[ 2 * counter_triangles + t ][2] = v22;
//       
//       data_triangle_vertices[ 3 * counter_triangles + t ][0] = v01;
//       data_triangle_vertices[ 3 * counter_triangles + t ][1] = v02;
//       data_triangle_vertices[ 3 * counter_triangles + t ][2] = v12;
//       
//     }
    
    
    // TODO: Check code until here.
    /**** Triangle -- Edge Relation ****/
    
    /* for each old edge, set the new first parent triangle */
    /* for each new edge, set the new first parent triangle */
    // checked
    for( int e = 0; e < counter_edges; e++ )
    {
      int p = data_edge_firstparent_triangle[e];
      
      assert( p != nullindex );
      
      int ei        = nullindex;
      int nfp_back  = nullindex;
      int nfp_front = nullindex;
      
      if( data_triangle_edges[p][0] == e ){
        ei = 0; nfp_back = 0; nfp_front = 1;
      } else if( data_triangle_edges[p][1] == e ) {
        ei = 1; nfp_back = 0; nfp_front = 2;
      } else if( data_triangle_edges[p][2] == e ) {
        ei = 2; nfp_back = 1; nfp_front = 2;
      } else 
        assert(false);
      
      assert( ei != nullindex );
      assert( data_triangle_edges[p][0] == e || data_triangle_edges[p][1] == e || data_triangle_edges[p][2] == e );
      assert( data_triangle_edges[p][ei] == e );
      
      data_edge_firstparent_triangle[ 0 * counter_edges + e ] = nfp_back  * counter_triangles + p;
      data_edge_firstparent_triangle[ 1 * counter_edges + e ] = nfp_front * counter_triangles + p;
      
    }
    
    
    /* for each triangle, relocate the data of the old edges' parent triangle */
    /* additionally, set the new parents */
    // checked
    for( int t  = 0; t  < counter_triangles;  t++ )
    for( int ei = 0; ei <                 3; ei++ )
    {
      int e = data_triangle_edges[t][ei];
      
      int t_back  = nullindex;
      int t_front = nullindex;
      int e_back  = nullindex;
      int e_front = nullindex;
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
      
      if( ei == 0 ){
        t_back = 0; t_front = 1; e_back = 0; e_front = 0; 
      } else if( ei == 1 ) {
        t_back = 0; t_front = 2; e_back = 1; e_front = 1; 
      } else if( ei == 2 ) {
        t_back = 1; t_front = 2; e_back = 2; e_front = 2; 
      } else 
        assert(false);
      
      
      int q = data_triangle_nextparents_of_edges[t][ei];
      
      if( q == nullindex ) {
        
        data_triangle_nextparents_of_edges[ t_back  * counter_triangles + t ][ e_back  ] = nullindex;
        data_triangle_nextparents_of_edges[ t_front * counter_triangles + t ][ e_front ] = nullindex;
        
      } else if( q != nullindex ) {
        
        int q_ei        = nullindex;
        int q_nfp_back  = nullindex;
        int q_nfp_front = nullindex;
        
        if(        data_triangle_edges[q][0] == e ){
          q_ei = 0; q_nfp_back = 0; q_nfp_front = 1;
        } else if( data_triangle_edges[q][1] == e ) {
          q_ei = 1; q_nfp_back = 0; q_nfp_front = 2;
        } else if( data_triangle_edges[q][2] == e ) {
          q_ei = 2; q_nfp_back = 1; q_nfp_front = 2;
        } else 
          assert(false);
        
        assert( q_ei != nullindex );
        assert( data_triangle_edges[q][0] == e || data_triangle_edges[q][1] == e || data_triangle_edges[q][2] == e );
        assert( data_triangle_edges[q][q_ei] == e );
        
        data_triangle_nextparents_of_edges[ t_back  * counter_triangles + t ][ e_back  ] = q_nfp_back  * counter_triangles + q;
        data_triangle_nextparents_of_edges[ t_front * counter_triangles + t ][ e_front ] = q_nfp_front * counter_triangles + q;
        
      } 
      
    }
    
    /* for each triangle, run over the new edges and add firstparents and parents */
    // checked 
    for( int t = 0; t < counter_triangles; t++ )
    {
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
      // new edges [ 01 02 ], [ 01 12 ], [ 02 12 ]
      
      data_edge_firstparent_triangle[ 2 * counter_edges + 0 * counter_triangles + t ] = 3 * counter_triangles + t;
      data_edge_firstparent_triangle[ 2 * counter_edges + 1 * counter_triangles + t ] = 3 * counter_triangles + t;
      data_edge_firstparent_triangle[ 2 * counter_edges + 2 * counter_triangles + t ] = 3 * counter_triangles + t;
      
      data_triangle_nextparents_of_edges[ 3 * counter_triangles + t ][0] = 0 * counter_triangles + t;
      data_triangle_nextparents_of_edges[ 3 * counter_triangles + t ][1] = 1 * counter_triangles + t;
      data_triangle_nextparents_of_edges[ 3 * counter_triangles + t ][2] = 2 * counter_triangles + t;
      
      data_triangle_nextparents_of_edges[ 0 * counter_triangles + t ][2] = nullindex;
      data_triangle_nextparents_of_edges[ 1 * counter_triangles + t ][1] = nullindex;
      data_triangle_nextparents_of_edges[ 2 * counter_triangles + t ][0] = nullindex;
      
    }
    
    
    
    /* for each new triangle, create the new vertices */
    // checked
    for( int t = 0; t < counter_triangles; t++ )
    {
      int v00 = data_triangle_vertices[t][0];
      int v11 = data_triangle_vertices[t][1];
      int v22 = data_triangle_vertices[t][2];
      
      int v01 = counter_vertices + data_triangle_edges[t][0];
      int v02 = counter_vertices + data_triangle_edges[t][1];
      int v12 = counter_vertices + data_triangle_edges[t][2];
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
      data_triangle_vertices[ 0 * counter_triangles + t ][0] = v00;
      data_triangle_vertices[ 0 * counter_triangles + t ][1] = v01;
      data_triangle_vertices[ 0 * counter_triangles + t ][2] = v02;
      
      data_triangle_vertices[ 1 * counter_triangles + t ][0] = v01;
      data_triangle_vertices[ 1 * counter_triangles + t ][1] = v11;
      data_triangle_vertices[ 1 * counter_triangles + t ][2] = v12;
      
      data_triangle_vertices[ 2 * counter_triangles + t ][0] = v02;
      data_triangle_vertices[ 2 * counter_triangles + t ][1] = v12;
      data_triangle_vertices[ 2 * counter_triangles + t ][2] = v22;
      
      data_triangle_vertices[ 3 * counter_triangles + t ][0] = v01;
      data_triangle_vertices[ 3 * counter_triangles + t ][1] = v02;
      data_triangle_vertices[ 3 * counter_triangles + t ][2] = v12;
      
    }
    
    /* for each new triangle, create the new edges */
    // checked
    for( int t = 0; t < counter_triangles; t++ )
    {
      int e00_01 = 0 * counter_edges + data_triangle_edges[t][0];
      int e00_02 = 0 * counter_edges + data_triangle_edges[t][1];
      int e01_11 = 1 * counter_edges + data_triangle_edges[t][0];
      int e02_22 = 1 * counter_edges + data_triangle_edges[t][1];
      int e11_12 = 0 * counter_edges + data_triangle_edges[t][2];
      int e12_22 = 1 * counter_edges + data_triangle_edges[t][2];
      
      int e01_02 = 2 * counter_edges + 0 * counter_triangles + t;
      int e01_12 = 2 * counter_edges + 1 * counter_triangles + t;
      int e02_12 = 2 * counter_edges + 2 * counter_triangles + t;
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
      data_triangle_edges[ 0 * counter_triangles + t ][0] = e00_01;
      data_triangle_edges[ 0 * counter_triangles + t ][1] = e00_02;
      data_triangle_edges[ 0 * counter_triangles + t ][2] = e01_02;
      
      data_triangle_edges[ 1 * counter_triangles + t ][0] = e01_11;
      data_triangle_edges[ 1 * counter_triangles + t ][1] = e01_12;
      data_triangle_edges[ 1 * counter_triangles + t ][2] = e11_12;
      
      data_triangle_edges[ 2 * counter_triangles + t ][0] = e02_12;
      data_triangle_edges[ 2 * counter_triangles + t ][1] = e02_22;
      data_triangle_edges[ 2 * counter_triangles + t ][2] = e12_22;
      
      data_triangle_edges[ 3 * counter_triangles + t ][0] = e01_02;
      data_triangle_edges[ 3 * counter_triangles + t ][1] = e01_12;
      data_triangle_edges[ 3 * counter_triangles + t ][2] = e02_12;
      
    }
    
    
    
    
    /* update the counters */
    
    counter_vertices  = counter_vertices + counter_edges;
    counter_edges     = 2 * counter_edges + 3 * counter_triangles;
    counter_triangles = 4 * counter_triangles;
    
    
    /* DONE */
    
    check();
}






FloatVector MeshSimplicial2D::get_triangle_midpoint( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = (   getcoordinates().getdata( get_triangle_vertices(t)[0], d ) 
                 + getcoordinates().getdata( get_triangle_vertices(t)[1], d ) 
                 + getcoordinates().getdata( get_triangle_vertices(t)[2], d )
                ) / 3.;
    return mid;
}

FloatVector MeshSimplicial2D::get_edge_midpoint    ( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = (   getcoordinates().getdata( get_edge_vertices(e)[0], d ) 
                 + getcoordinates().getdata( get_edge_vertices(e)[1], d )
                ) / 2.;
    return mid;
}




//         if( data_edge_firstparent_triangle[ t_e2 ] == t_old ) {
//           
//           data_edge_firstparent_triangle[ t_e2 ] = counter_triangles + ot;
//           
//         } else {
//           
//           int current_edge = data_vertex_firstparent_edge[ e_front_vertex ];
//           while( data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] != e )
//             current_edge = data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ];
//           data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] = counter_edge;
//           
//         }
//       
//         /* TODO: triangle parents of opposing vertex */
//         if( data_vertex_firstparent_triangle[ t_v2 ] == t_old ) { 
//           int nextparent_tri = data_triangle_nextparents_of_vertices[ current_tri ][ 2 ];
//           data_vertex_firstparent_triangle[ t_v2 ] = counter_edges;
//         } else {
//           int current_tri = data_vertex_firstparent_triangle[ t_v2 ];
//           while( data_triangle_nextparents_of_vertices[ current_tri ][ 2 ] != t_old )
//             current_tri = data_triangle_nextparents_of_vertices[ current_tri ][ 2 ];
//           data_triangle_nextparents_of_vertices[ current_tri ][ 2 ] = counter_edge;
//         }
