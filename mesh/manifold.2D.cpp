
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
#include "manifold.2D.hpp"




/* from local edge index, get the local vertex indices */

const std::array<int,2> localedge_to_localvertices( int ei )
{
  std::array<int,2> ret;
  if( ei == 0 ) {
    ret[0] = 0; ret[1] = 1;
  } else
  if( ei == 1 ) {
    ret[0] = 0; ret[1] = 2;
  } else
  if( ei == 2 ) {
    ret[0] = 1; ret[1] = 2;
  } else
    assert( false );
  return ret;
}












ManifoldTriangulation2D::ManifoldTriangulation2D( int outerdim )
:
    outerdimension(outerdim),
    coordinates(outerdim,0),
    
    counter_triangles(0),
    counter_edges(0),
    counter_vertices(0),
    
    data_triangle_vertices(0),
    data_triangle_edges(0),
    data_edge_parents(0)
//     data_vertex_firstparent(0)   
{
    check();
}


ManifoldTriangulation2D::ManifoldTriangulation2D( 
    int outerdim,
    const Coordinates& coords,
    const std::vector<std::array<int,3>> triangle_vertex_list
)
:
    outerdimension(outerdim),
    coordinates(outerdim,0),
    
    counter_triangles(0),
    counter_edges(0),
    counter_vertices(0),
    
    data_triangle_vertices(0),
    data_triangle_edges(0),
    data_edge_parents(0)
//     data_vertex_firstparent(0)
{
    check();
}


ManifoldTriangulation2D::~ManifoldTriangulation2D()
{
    
}


void ManifoldTriangulation2D::check() const
{
    // FIXME: add code 
    /* checke welche listen vorhanden sein sollen. */
    /* check die alle subsimplex listen. Dabei auch die supersimplex liste prüfen */
    /* Supersimplex listen checken */
}

void ManifoldTriangulation2D::print( std::ostream& os ) const
{
    os << "Printe Triangulation of 2D Manifold!" << std::endl;
}







/* Count number of elements */

int ManifoldTriangulation2D::getouterdimension() const
{
    return outerdimension;
}

Coordinates& ManifoldTriangulation2D::getcoordinates()
{
    return coordinates;
}


const Coordinates& ManifoldTriangulation2D::getcoordinates() const
{
    return coordinates;
}







int ManifoldTriangulation2D::count_triangles() const
{
    return counter_triangles;
}

int ManifoldTriangulation2D::count_edges() const
{
    return counter_edges;
}

int ManifoldTriangulation2D::count_vertices() const
{
    return counter_vertices;
}







bool ManifoldTriangulation2D::contains_triangle_edge  ( int t, int e ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= e && e < counter_edges );
    return    ( data_triangle_edges[t][0] == e )
           || ( data_triangle_edges[t][1] == e )
           || ( data_triangle_edges[t][2] == e )
           ;
} 

bool ManifoldTriangulation2D::contains_triangle_vertex( int t, int v ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= v && v < counter_vertices );
    return    ( data_triangle_vertices[t][0] == v )
           || ( data_triangle_vertices[t][1] == v )
           || ( data_triangle_vertices[t][2] == v )
           ;
} 

bool ManifoldTriangulation2D::contains_edge_vertex    ( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    
    int t = get_edge_parents( e )[0];
    assert( contains_triangle_edge( t, e ) );
    
    int ei = indexof_triangle_edge( t, e );
    std::array<int,2> vi = localedge_to_localvertices( ei );
    
    return    ( data_triangle_vertices[t][vi[0]] == v )
           || ( data_triangle_vertices[t][vi[1]] == v )
           ;
} 







int ManifoldTriangulation2D::indexof_triangle_edge   ( int t, int e ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= e && e < counter_edges );
    if ( data_triangle_edges[t][0] == e )
      return 0;
    else
    if ( data_triangle_edges[t][1] == e )
      return 1;
    else
    if ( data_triangle_edges[t][2] == e )
      return 2;
    else
      assert(false);
} 

 
int ManifoldTriangulation2D::indexof_triangle_vertex( int t, int v ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= v && v < counter_vertices );
    if ( data_triangle_vertices[t][0] == v )
      return 0;
    else
    if ( data_triangle_vertices[t][1] == v )
      return 1;
    else
    if ( data_triangle_vertices[t][2] == v )
      return 2;
    else
      assert(false);
} 


int ManifoldTriangulation2D::indexof_edge_vertex    ( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    
    int t = get_edge_parents( e )[0];
    assert( contains_triangle_edge( t, e ) );
    
    int ei = indexof_triangle_edge( t, e );
    std::array<int,2> vi = localedge_to_localvertices( ei );
    
    if ( data_triangle_vertices[t][vi[0]] == v )
      return vi[0];
    else
    if ( data_triangle_vertices[t][vi[1]] == v )
      return vi[1];
    else
      assert(false); 
} 










const std::array<int,3> ManifoldTriangulation2D::get_triangle_edges   ( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    return data_triangle_vertices[t];
} 

const std::array<int,3> ManifoldTriangulation2D::get_triangle_vertices( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    return data_triangle_edges[t];
}

const std::array<int,2> ManifoldTriangulation2D::get_edge_vertices    ( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    
    int t = get_edge_parents( e )[0];
    assert( contains_triangle_edge( t, e ) );
    
    std::array<int,2> vi = localedge_to_localvertices( e );
    std::array<int,2> vs;
    vs[0] = data_triangle_vertices[t][vi[0]];
    vs[1] = data_triangle_vertices[t][vi[0]];
    
    return vs;
}








/* is a neighbor? */
bool ManifoldTriangulation2D::is_triangle_neighbor( int t, int nt ) const
{
    assert( 0 <=  t &&  t < counter_triangles );
    assert( 0 <= nt && nt < counter_triangles );
    
    if( t == nt ) return false;
    
    for( int ei = 0; ei < 3; ei++ ) 
    for( int pi = 0; pi < 2; pi++ ) 
    {
        int e = data_triangle_edges[t][ei];
        int p = data_edge_parents[e][pi]; 
        if( p == nt ) return true;
    }
    return false;
}

/* get index of a neighbor */
int ManifoldTriangulation2D::indexof_triangle_neighbor( int t, int nt ) const
{
    assert( 0 <=  t &&  t < counter_triangles );
    assert( 0 <= nt && nt < counter_triangles );
    
    if( t == nt ) return false;
    
    for( int ei = 0; ei < 3; ei++ ) 
    for( int pi = 0; pi < 2; pi++ ) 
    {
        int e = data_triangle_edges[t][ei];
        int p = data_edge_parents[e][pi]; 
        if( p == nt ) return pi;
    }
    return false;
} 

/* get the triangle neighbors */
const std::array<int,3> ManifoldTriangulation2D::get_triangle_neighbors( int t ) const 
{
    assert( 0 <=  t &&  t < counter_triangles );
    
    std::array<int,3> ret;
    
    for( int ei = 0; ei < 3; ei++ ) 
    {
        int e = data_triangle_edges[t][ei];
        
        if( data_edge_parents[e][0] != t )
          ret[ei] = data_edge_parents[e][1];
        else if( data_edge_parents[e][1] != t )
          ret[ei] = data_edge_parents[e][0];
        else 
          assert(false);
      
    }
    return ret;
} 





bool ManifoldTriangulation2D::is_edge_between( int t1, int t2, int e ) const
{
    assert( 0 <= t1 && t1 < counter_triangles );
    assert( 0 <= t2 && t2 < counter_triangles );
    
    assert( 0 <= e && e < counter_edges );
    
    if( ! is_triangle_neighbor( t1, t2 ) ) return false;
    if( ! contains_triangle_edge( t1, e ) ) return false;
    if( ! contains_triangle_edge( t2, e ) ) return false;
    
}

int ManifoldTriangulation2D::get_edge_between( int t1, int t2 ) const
{
    assert( 0 <= t1 && t1 < counter_triangles );
    assert( 0 <= t2 && t2 < counter_triangles );
    
    int e;
    
    
    assert( is_edge_between( t1, t2, e ) );
    
    return e;
}


















/* parents of an edge */

bool ManifoldTriangulation2D::is_edge_parent( int t, int e ) const
{
    assert( 0 <=  t &&  t < counter_triangles );
    assert( 0 <=  e &&  e < counter_edges );
    
    return ( data_edge_parents[e][0] == t ) || ( data_edge_parents[e][1] == t );
}


int ManifoldTriangulation2D::indexof_edge_parent( int t, int e ) const
{
    assert( 0 <=  t &&  t < counter_triangles );
    assert( 0 <=  e &&  e < counter_edges );
    
    if( data_edge_parents[e][0] == t )
      return 0;
    else
    if( data_edge_parents[e][1] == t )
      return 1;
    else
      assert(false);
}


int ManifoldTriangulation2D::count_edge_parents( int e ) const 
{
    assert( 0 <=  e &&  e < counter_edges );
    
    if( data_edge_parents[e][0] != nullindex && data_edge_parents[e][1] != nullindex )
      return 2;
    else
    if( data_edge_parents[e][0] != nullindex && data_edge_parents[e][1] == nullindex )
      return 1;
    else
    if( data_edge_parents[e][0] == nullindex && data_edge_parents[e][1] != nullindex )
      return 1;
    else 
      assert(false);
    
}

const std::array<int,2> ManifoldTriangulation2D::get_edge_parents( int e ) const
{
    assert( 0 <=  e &&  e < counter_edges );
    
    return data_edge_parents[e];
    
}











// // int ManifoldTriangulation2D::is_vertexparents_cyclic( int v );
// // 
// // int ManifoldTriangulation2D::count_vertexparents( int v, int* ts );
// // int ManifoldTriangulation2D::list_vertexparents( int v, int* ts );
// // int ManifoldTriangulation2D::count_vertexparents_edges( int v, int* es );
// // int ManifoldTriangulation2D::list_vertexparents_edges( int v, int* es );



/* bisect a single edge */

void ManifoldTriangulation2D::bisect_edge( int e )
{
  
}



