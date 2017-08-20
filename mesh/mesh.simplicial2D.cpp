
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
        
        assert( data_edge_vertices[t1][0] != data_edge_vertices[t2][0] || data_edge_vertices[t1][1] != data_edge_vertices[t2][1] || data_edge_vertices[t1][2] != data_edge_vertices[t2][2] );
        assert( data_edge_vertices[t1][0] != data_edge_vertices[t2][0] || data_edge_vertices[t1][1] != data_edge_vertices[t2][2] || data_edge_vertices[t1][2] != data_edge_vertices[t2][1] );
        assert( data_edge_vertices[t1][0] != data_edge_vertices[t2][1] || data_edge_vertices[t1][1] != data_edge_vertices[t2][0] || data_edge_vertices[t1][2] != data_edge_vertices[t2][2] );
        assert( data_edge_vertices[t1][0] != data_edge_vertices[t2][1] || data_edge_vertices[t1][1] != data_edge_vertices[t2][2] || data_edge_vertices[t1][2] != data_edge_vertices[t2][0] );
        assert( data_edge_vertices[t1][0] != data_edge_vertices[t2][2] || data_edge_vertices[t1][1] != data_edge_vertices[t2][0] || data_edge_vertices[t1][2] != data_edge_vertices[t2][1] );
        assert( data_edge_vertices[t1][0] != data_edge_vertices[t2][2] || data_edge_vertices[t1][1] != data_edge_vertices[t2][1] || data_edge_vertices[t1][2] != data_edge_vertices[t2][0] );
        
        
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
    
    for( const auto& triple : data_triangle_vertices )
      std::cout << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    os << "Edge first parent triangles" << std::endl;
    
    for( int fp : data_edge_firstparent_triangle )
      std::cout << fp << nl;
    
    os << "Triangle next parents of edges" << std::endl;
    
    for( const auto& triple : data_edge_nextparents_of_vertices )
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






bool MeshSimplicial2D::dimensioncounted( int dim ) const
{
    assert( 0 <= dim && dim <= 2 );
    return true;
}

int MeshSimplicial2D::countsimplices( int dim ) const
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
    assert( 0 <= sub && sub < sup && sup <= 2 );
    return true;
}

IndexMap MeshSimplicial2D::getsubsimplices( int sup, int sub, int cell ) const
{
  
  if( sup == 2 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_triangles() );
    auto temp = get_triangle_edges(cell);
    return IndexMap( IndexRange(0,2), IndexRange( 0, count_edges()-1 ), std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 2 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_triangles() );
    auto temp = get_triangle_vertices(cell);
    return IndexMap( IndexRange(0,2), IndexRange( 0, count_vertices()-1 ) , std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 1 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    auto temp = get_edge_vertices(cell);
    return IndexMap( IndexRange(0,1), IndexRange( 0, count_vertices()-1 ) , std::vector<int>( temp.begin(), temp.end() ) );
    
  } else {
    
    assert(false);
    
  }
   
}

bool MeshSimplicial2D::supersimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub < sup && sup <= 2 );
    return true;
}

const std::vector<int> MeshSimplicial2D::getsupersimplices( int sup, int sub, int cell ) const
{
  
  assert( 1 == sup );
  assert( 0 == sub );
  assert( 0 <= cell && cell < count_vertices() );
  
  if( sup == 2 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    auto temp = get_triangle_parents_of_edge( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 2 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    auto temp = get_triangle_parents_of_vertex( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 1 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    auto temp = get_edge_parents_of_vertex( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
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





void MeshSimplicial2D::bisect_edge( int e )
{
    // TODO: Rewrite this section for the 2D case 
    
    assert( 0 <= e && e < counter_edges );
    check();
    
    /* Collect the old data */
    
    int vertex_back  = data_edge_vertices[e][0];
    int vertex_front = data_edge_vertices[e][1];
    int nextparents_back  = data_edge_nextparents_of_vertices[e][0];
    int nextparents_front = data_edge_nextparents_of_vertices[e][1];
    int firstparent_back  = data_vertex_firstparent_edge[vertex_back ];
    int firstparent_front = data_vertex_firstparent_edge[vertex_front];
    
    int back_previousparent   = nullindex;
    int back_previousparent_localindex = nullindex;
    if( e != get_vertex_firstparent_edge( vertex_back ) )
      for( back_previousparent = get_vertex_firstparent_edge( vertex_back  );
           back_previousparent != nullindex && get_vertex_nextparent_edge(vertex_back,back_previousparent) != e; 
           back_previousparent = get_vertex_nextparent_edge(vertex_back,back_previousparent) 
         ); 
    if( back_previousparent  != nullindex ) 
      back_previousparent_localindex  = indexof_edge_vertex( back_previousparent,  vertex_back  );
    
    int front_previousparent  = nullindex;
    int front_previousparent_localindex = nullindex;
    if( e != get_vertex_firstparent_edge( vertex_front ) )
      for( front_previousparent = get_vertex_firstparent_edge( vertex_front );
           front_previousparent != nullindex && get_vertex_nextparent_edge(vertex_front,front_previousparent) != e; 
           front_previousparent = get_vertex_nextparent_edge(vertex_front,front_previousparent) 
         ); 
    if( front_previousparent != nullindex )
      front_previousparent_localindex = indexof_edge_vertex( front_previousparent, vertex_front );
    
    /* Assemble the data */
    
    FloatVector midcoordinate = get_edge_midpoint( e );
    
    int ne = counter_edges;
    int nv = counter_vertices;
    
    int back_backvertex       = vertex_back;
    int back_frontvertex      = nv;
    int front_backvertex      = nv;
    int front_frontvertex     = vertex_front;
    
    int back_backnextparent   = nextparents_back;
    int back_frontnextparent  = ne;
    int front_backnextparent  = nullindex;
    int front_frontnextparent = nextparents_front;
    
    int firstparent_newvertex = nv;
    
    /* Allocate memory */
    
    data_edge_nextparents_of_vertices.resize  ( counter_edges    + 1 );
    data_edge_vertices.resize     ( counter_edges    + 1 );
    data_vertex_firstparent_edge.resize( counter_vertices + 1 );
    
    /* Write in the data */
    
    data_edge_vertices[e ][0] = back_backvertex;
    data_edge_vertices[e ][1] = back_frontvertex;
    data_edge_vertices[ne][0] = front_backvertex;
    data_edge_vertices[ne][1] = front_frontvertex;
    
    data_edge_nextparents_of_vertices[e ][0] = back_backnextparent;
    data_edge_nextparents_of_vertices[e ][1] = back_frontnextparent;
    data_edge_nextparents_of_vertices[ne][0] = front_backnextparent;
    data_edge_nextparents_of_vertices[ne][1] = front_frontnextparent;
    
    data_vertex_firstparent_edge[nv] = e;
    
    if( back_previousparent  != nullindex ) {
      assert( data_edge_nextparents_of_vertices[ back_previousparent ][ back_previousparent_localindex ] == e );
      data_edge_nextparents_of_vertices[ back_previousparent ][ back_previousparent_localindex ] = e;
    } else {
      assert( data_vertex_firstparent_edge[ vertex_back ] == e );
      data_vertex_firstparent_edge[ vertex_back ] = e;
    }
    
    if( front_previousparent != nullindex ) {
      assert( data_edge_nextparents_of_vertices[ front_previousparent ][ front_previousparent_localindex ] == ne );
      data_edge_nextparents_of_vertices[ front_previousparent ][ front_previousparent_localindex ] = ne;
    } else {
      assert( data_vertex_firstparent_edge[ vertex_front ] == e );
      data_vertex_firstparent_edge[ vertex_front ] = ne;
    }    
    
    getcoordinates().append( midcoordinate );
        
    /* Update counter */
    counter_edges++;
    counter_vertices++;
    
    /* Done */
    
    
    check();
    
}


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




