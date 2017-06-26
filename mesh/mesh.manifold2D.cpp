
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
#include "mesh.manifold2D.hpp"




MeshManifold2D::MeshManifold2D( int outerdim )
:
    Mesh( 2, outerdim ),
    
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


MeshManifold2D::MeshManifold2D( 
    int outerdim,
    const Coordinates& coords,
    const std::vector<std::array<int,3>> triangle_vertex_list
)
:
    Mesh( 2, outerdim ),
    
    counter_triangles( triangle_vertex_list.size() ),
    counter_edges(0),
    counter_vertices(0),
    
    data_triangle_vertices( triangle_vertex_list ),
    data_triangle_edges( triangle_vertex_list.size(), {nullindex,nullindex,nullindex} ),
    data_edge_parents(0)
//     data_vertex_firstparent(0)
{
    
    getcoordinates() = coords;
    
    /* 1. Count triangles, transfer data */ 
    /* DONE */
    
    
    /* 2. Count vertices */
    counter_vertices = 0;
    for( const auto& tuple : data_triangle_vertices )
    for( const int& vertex : tuple )
      counter_vertices = counter_vertices < vertex ? vertex : counter_vertices; 
//       maximum( counter_vertices, vertex );
    counter_vertices += 1;
      
    /* 3. Count edges, T->E, E->T  */
    /* FIXME: quadratic run-time. Can be pushed down to linearlogarithmic. */
    
    /* first the inner edges */
    for( int t1 = 0; t1 < counter_triangles; t1++ )
    for( int t2 = 0; t2 < t1               ; t2++ )
    for( int e1 = 0; e1 < 3; e1++ )
    for( int e2 = 0; e2 < 3; e2++ )
    {
        const auto& ve1 = duple_from_triple( data_triangle_vertices[t1], edgeindex_to_vertexindices( e1 ) );
        const auto& ve2 = duple_from_triple( data_triangle_vertices[t2], edgeindex_to_vertexindices( e2 ) );
        if( vertexpairs_equivalent( ve1, ve2 ) )
        {
            data_edge_parents.push_back( { t1, t2 } );
            data_triangle_edges[t1][ e1 ] = data_edge_parents.size() - 1;
            data_triangle_edges[t2][ e2 ] = data_edge_parents.size() - 1;
        }
    }
    
    /* then the outer edges */
    for( int t = 0; t < counter_triangles; t++ )
    for( int e = 0; e < 3; e++ )
    {
        if( data_triangle_edges[t][e] == nullindex ) 
        {
            data_edge_parents.push_back( { t, nullindex } );
            data_triangle_edges[t][e] = data_edge_parents.size() - 1;
        }
    }
    
    counter_edges = data_edge_parents.size();
    
    check();
}


MeshManifold2D::MeshManifold2D( 
    int outerdim,
    const Coordinates& coords,
    const std::vector<std::array<int,3>> triangle_vertices,
    const std::vector<std::array<int,3>> triangle_edges,
    const std::vector<std::array<int,2>> edge_parents
)
:
    Mesh( 2, outerdim ),
    
    counter_triangles( triangle_vertices.size() ),
    counter_edges( edge_parents.size() ),
    counter_vertices( coords.getnumber() ),
    
    data_triangle_vertices( triangle_vertices ),
    data_triangle_edges( triangle_edges ),
    data_edge_parents( edge_parents )
//     data_vertex_firstparent(0)
{
    
    getcoordinates() = coords;
    
    check();
}


MeshManifold2D::~MeshManifold2D()
{
    
}


void MeshManifold2D::check() const
{
    /* * Data integrity
     * 
     * - each triangle: each vertex is non-null
     * - each triangle: each vertex is unique
     * - each triangle: each edge is non-null
     * - each triangle: each edge is unique
     * 
     */
    
//     cout << count_vertices() << space << getcoordinates().getnumber() << nl;

    assert( counter_triangles == data_triangle_vertices.size() );
    assert( counter_triangles == data_triangle_edges.size() );
    assert( counter_edges     == data_edge_parents.size() );
    assert( count_vertices() == getcoordinates().getnumber() );
    
    for( int t = 0; t < counter_triangles; t++ )
    {
        assert( data_triangle_vertices[t][0] != nullindex );
        assert( data_triangle_vertices[t][1] != nullindex );
        assert( data_triangle_vertices[t][2] != nullindex );
        
        assert( data_triangle_vertices[t][0] != data_triangle_vertices[t][1] );
        assert( data_triangle_vertices[t][1] != data_triangle_vertices[t][2] );
        assert( data_triangle_vertices[t][2] != data_triangle_vertices[t][0] );
        
        assert( data_triangle_edges[t][0] != nullindex );
        assert( data_triangle_edges[t][1] != nullindex );
        assert( data_triangle_edges[t][2] != nullindex );
        
        assert( data_triangle_edges[t][0] != data_triangle_edges[t][1] );
        assert( data_triangle_edges[t][1] != data_triangle_edges[t][2] );
        assert( data_triangle_edges[t][2] != data_triangle_edges[t][0] );
        
    }
    
    /* * Data integrity
     * 
     * - each edge: the parents are unique
     * - each edge: at least one parent is non-null 
     * - each edge: the first parent is non-null
     * - each edge: each parent triangle points to the edge
     * 
     */
    
    assert( counter_edges == data_edge_parents.size() );
    
    for( int e = 0; e < counter_edges; e++ )
    {
        
//         std::cout << "vergleich " << e << " : " << data_edge_parents[e][0] << space << data_edge_parents[e][1] << nl;
        
        assert( data_edge_parents[e][0] != data_edge_parents[e][1] );
        assert( data_edge_parents[e][0] != nullindex || data_edge_parents[e][1] != nullindex );
        assert( data_edge_parents[e][0] != nullindex );
        
        int p0 = data_edge_parents[e][0];
        
        assert( data_triangle_edges[p0][0] == e || data_triangle_edges[p0][1] == e || data_triangle_edges[p0][2] == e );
        
        {
          std::array<int,2> e_vertices1 = get_edge_vertices( e );
          std::array<int,2> e_vertices2 = duple_from_triple( data_triangle_vertices[p0], edgeindex_to_vertexindices( indexof_triangle_edge( p0, e ) ) );
          
          assert( vertexpairs_equivalent( e_vertices1, e_vertices2 ) );
        }
        
        int p1 = data_edge_parents[e][1];
        
        if( p1 == nullindex) continue; 
        
        assert( data_triangle_edges[p1][0] == e || data_triangle_edges[p1][1] == e || data_triangle_edges[p1][2] == e );
        
        {
          std::array<int,2> e_vertices1 = get_edge_vertices( e );
          std::array<int,2> e_vertices2 = duple_from_triple( data_triangle_vertices[p1], edgeindex_to_vertexindices( indexof_triangle_edge( p1, e ) ) );
          
          assert( vertexpairs_equivalent( e_vertices1, e_vertices2 ) );
        }
        
    }
    
    /* check for duplicates */
    
    for( int t1 = 0; t1 < counter_triangles; t1++ )
    for( int t2 = 0; t2 < t1; t2++ )
    {
      assert( ! triple_equivalent( data_triangle_edges[t1], data_triangle_edges[t2] ) );
      assert( ! triple_equivalent( data_triangle_vertices[t1], data_triangle_vertices[t2] ) );
    }
    
//     for( int e1 = 0; e1 < counter_edges; e1++ )
//     for( int e2 = 0; e2 < e1; e2++ )
//       if     ( data_edge_parents[e1][0] == nullindex && data_edge_parents[e2][0] == nullindex )
//         assert( data_edge_parents[e1][1] != data_edge_parents[e2][1] );
//       else if( data_edge_parents[e1][0] == nullindex && data_edge_parents[e2][1] == nullindex )
//         assert( data_edge_parents[e1][1] != data_edge_parents[e2][0] );
//       else if( data_edge_parents[e1][1] == nullindex && data_edge_parents[e2][0] == nullindex )
//         assert( data_edge_parents[e1][0] != data_edge_parents[e2][1] );
//       else if( data_edge_parents[e1][1] == nullindex && data_edge_parents[e2][1] == nullindex )
//         assert( data_edge_parents[e1][0] != data_edge_parents[e2][0] );
//       else
//         assert( ! duple_equivalent( data_edge_parents[e1], data_edge_parents[e2] ) );
    
}






void MeshManifold2D::print( std::ostream& os ) const
{
    os << "Printe Triangulation of 2D Manifold!" << std::endl;
    
    os << counter_triangles << space << counter_edges << space << counter_vertices << nl;
    
    os << "Triangle vertices" << std::endl;
    
    for( const auto& tuple : data_triangle_vertices )
      std::cout << tuple[0] << space << tuple[1] << space << tuple[2] << nl;
    
    os << "Triangle edges" << std::endl;
    
    for( const auto& tuple : data_triangle_edges )
      std::cout << tuple[0] << space << tuple[1] << space << tuple[2] << nl;
    
    os << "Edge parents " << std::endl;
    
    for( const auto& tuple : data_edge_parents )
      std::cout << tuple[0] << space << tuple[1] << nl;
    
    os << "Edge vertices" << std::endl;
    
    for( int e = 0; e < counter_edges; e++ )
      std::cout << e << ":" << get_edge_vertices(e)[0] << space << get_edge_vertices(e)[1] << nl;
    
    os << "Finished printing" << nl;
    
}






bool MeshManifold2D::dimensioncounted( int ) const
{
    return true;
}

int MeshManifold2D::countsimplices( int dim ) const
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

bool MeshManifold2D::subsimplices_listed( int, int ) const
{
  return true;
}

const IndexMap MeshManifold2D::getsubsimplices( int sup, int sub, int cell ) const
{
  assert( 0 <= sup && sup <= 2 );
  assert( 0 <= sub && sub <= 2 );
  assert( sup != sub );
  assert( sub < sup );
  if( sup == 2 ) assert( 0 <= cell && cell < count_triangles() );
  if( sup == 1 ) assert( 0 <= cell && cell < count_edges()     );
  if( sup == 0 ) assert( 0 <= cell && cell < count_vertices()  );
  
  if ( sup == 2 && sub == 1 ) {
    
    auto temp = get_triangle_edges(cell);
    return IndexMap( IndexRange(0,2), std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 2 && sub == 0 ) {
    
    auto temp = get_triangle_vertices(cell);
    return IndexMap( IndexRange(0,2), std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 1 && sub == 0 ) {
    
    auto temp = get_edge_vertices(cell);
    return IndexMap( IndexRange(0,1), std::vector<int>( temp.begin(), temp.end() ) );
    
  } else {
    
    assert(false);
    
  }
    
}

bool MeshManifold2D::supersimplices_listed( int, int ) const
{
  return true;
}

const std::vector<int> MeshManifold2D::getsupersimplices( int sup, int sub, int cell ) const
{
  assert( 0 <= sup && sup <= 2 );
  assert( 0 <= sub && sub <= 2 );
  assert( sup != sub );
  assert( sub < sup );
  if( sub == 2 ) assert( 0 <= cell && cell < count_triangles() );
  if( sub == 1 ) assert( 0 <= cell && cell < count_edges()     );
  if( sub == 0 ) assert( 0 <= cell && cell < count_vertices()  );
  
  if ( sup == 2 && sub == 1 ) {
    
    auto temp = get_triangle_edge_parents( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 2 && sub == 0 ) {
    
    auto temp = get_triangle_parents_of_vertex( cell );
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 1 && sub == 0 ) {
    
    auto temp = get_edge_parents_of_vertex( cell ); 
    return std::vector<int>( temp.begin(), temp.end() );
    
  } else {
    
    assert(false);
    
  }
  
}








/* Count number of elements */


int MeshManifold2D::count_triangles() const
{
    assert( counter_triangles == data_triangle_edges.size()    );
    assert( counter_triangles == data_triangle_vertices.size() );
    return counter_triangles;
}

int MeshManifold2D::count_edges() const
{
    return counter_edges;
}

int MeshManifold2D::count_vertices() const
{
    return counter_vertices;
}







bool MeshManifold2D::contains_triangle_edge  ( int t, int e ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= e && e < counter_edges );
    return    ( data_triangle_edges[t][0] == e )
           || ( data_triangle_edges[t][1] == e )
           || ( data_triangle_edges[t][2] == e )
           ;
} 

bool MeshManifold2D::contains_triangle_vertex( int t, int v ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= v && v < counter_vertices );
    return    ( data_triangle_vertices[t][0] == v )
           || ( data_triangle_vertices[t][1] == v )
           || ( data_triangle_vertices[t][2] == v )
           ;
} 

bool MeshManifold2D::contains_edge_vertex    ( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    
    int t = get_triangle_edge_parents( e )[0];
    assert( is_not_nullindex( t ) );
    assert( contains_triangle_edge( t, e ) );
    
    int ei = indexof_triangle_edge( t, e );
    std::array<int,2> vi = edgeindex_to_vertexindices( ei );
    
    return    ( data_triangle_vertices[t][vi[0]] == v )
           || ( data_triangle_vertices[t][vi[1]] == v )
           ;
} 







int MeshManifold2D::indexof_triangle_edge   ( int t, int e ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= e && e < counter_edges );
    if     ( data_triangle_edges[t][0] == e ) return 0;
    else if( data_triangle_edges[t][1] == e ) return 1;
    else if( data_triangle_edges[t][2] == e ) return 2;
    else                                      assert(false);
} 

int MeshManifold2D::indexof_triangle_vertex( int t, int v ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= v && v < counter_vertices );
    if     ( data_triangle_vertices[t][0] == v ) return 0;
    else if( data_triangle_vertices[t][1] == v ) return 1;
    else if( data_triangle_vertices[t][2] == v ) return 2;
    else                                         assert(false);
} 

int MeshManifold2D::indexof_edge_vertex    ( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    
    int t = get_triangle_edge_parents( e )[0];
    assert( contains_triangle_edge( t, e ) );
    
    int ei = indexof_triangle_edge( t, e );
    std::array<int,2> vi = edgeindex_to_vertexindices( ei );
    
    if     ( data_triangle_vertices[t][vi[0]] == v ) return vi[0];
    else if( data_triangle_vertices[t][vi[1]] == v ) return vi[1];
    else                                             assert(false); 
} 










const std::array<int,3> MeshManifold2D::get_triangle_edges   ( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    return data_triangle_edges[t];
} 

const std::array<int,3> MeshManifold2D::get_triangle_vertices( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    return data_triangle_vertices[t];
}

const std::array<int,2> MeshManifold2D::get_edge_vertices    ( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    
    int t = get_triangle_edge_parents( e )[0];
    assert( contains_triangle_edge( t, e ) );
    
    int el = indexof_triangle_edge( t, e );
    
    std::array<int,2> vi = edgeindex_to_vertexindices( el );
    
    return { data_triangle_vertices[t][ vi[0] ], data_triangle_vertices[t][ vi[1] ] };
}








/* is a neighbor? */
bool MeshManifold2D::is_triangle_neighbor( int t, int nt ) const
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
int MeshManifold2D::indexof_triangle_neighbor( int t, int nt ) const
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
const std::array<int,3> MeshManifold2D::get_triangle_neighbors( int t ) const 
{
    assert( 0 <=  t &&  t < counter_triangles );
    
    std::array<int,3> ret;
    
    for( int ei = 0; ei < 3; ei++ ) 
    {
        int e = data_triangle_edges[t][ei];
        
        if( data_edge_parents[e][0] == t )
          ret[ei] = data_edge_parents[e][1];
        else if( data_edge_parents[e][1] == t )
          ret[ei] = data_edge_parents[e][0];
        else 
          assert(false);
      
    }
    return ret;
} 





bool MeshManifold2D::is_edge_between( int t1, int t2, int e ) const
{
    assert( 0 <= t1 && t1 < counter_triangles );
    assert( 0 <= t2 && t2 < counter_triangles );
    
    assert( 0 <= e && e < counter_edges );
    
    if( ! is_triangle_neighbor( t1, t2 )  ) return false;
    if( ! contains_triangle_edge( t1, e ) ) return false;
    if( ! contains_triangle_edge( t2, e ) ) return false;
    
    return true;
}

int MeshManifold2D::get_edge_between( int t1, int t2 ) const
{
    assert( 0 <= t1 && t1 < counter_triangles );
    assert( 0 <= t2 && t2 < counter_triangles );
    
    assert( is_triangle_neighbor( t1, t2 ) );
    
    int n = indexof_triangle_neighbor( t1, t2 );
    
    return neighborindex_to_edgeindex( n );
}








int MeshManifold2D::get_forwarding_edge ( int t, int el, int o )
{
    assert( 0 <=  t &&  t < counter_triangles );
    assert( 0 <= el && el < 3 );
    assert( o == 1 || o == -1 );
    
    int e = data_triangle_edges[t][el];
    
    std::array<int,2> e_vertices = get_edge_vertices(e);
    
    int v = ( o == 1 ) ? e_vertices[0] : e_vertices[1];
    
    assert( contains_triangle_vertex( t, v ) );
    
    int v_local = indexof_triangle_vertex( t, v );
    
    return vertexindex_to_opposing_edgeindex( v_local );
}

int MeshManifold2D::get_backwarding_edge( int t, int el, int o )
{
    assert( 0 <=  t &&  t < counter_triangles );
    assert( 0 <= el && el < 3 );
    assert( o == 1 || o == -1 );
    
    int e = data_triangle_edges[t][el];
    
    std::array<int,2> e_vertices = get_edge_vertices(e);
    
    int v = ( o == 1 ) ? e_vertices[1] : e_vertices[0];
    
    assert( contains_triangle_vertex( t, v ) );
    
    int v_local = indexof_triangle_vertex( t, v );
    
    return vertexindex_to_opposing_edgeindex( v_local );
}



int MeshManifold2D::get_prev_edge    ( int el, int o )
{
    assert( 0 <= el && el < 3 );
    assert( o == 1 || o == -1 );
    if( el == 0 && o ==  1 ) return 2;
    if( el == 0 && o == -1 ) return 1;
    if( el == 1 && o ==  1 ) return 0;
    if( el == 1 && o == -1 ) return 2;
    if( el == 2 && o ==  1 ) return 1;
    if( el == 2 && o == -1 ) return 0;
    assert(false);
}

int MeshManifold2D::get_next_edge    ( int el, int o )
{
    assert( 0 <= el && el < 3 );
    assert( o == 1 || o == -1 );
    if( el == 0 && o ==  1 ) return 1;
    if( el == 0 && o == -1 ) return 2;
    if( el == 1 && o ==  1 ) return 2;
    if( el == 1 && o == -1 ) return 0;
    if( el == 2 && o ==  1 ) return 0;
    if( el == 2 && o == -1 ) return 1;
    assert(false);
}


// 01 02 12 -> 2 1 0

int MeshManifold2D::get_prev_neighbor( int nl, int o )
{
    assert( 0 <= nl && nl < 3 );
    assert( o == 1 || o == -1 );
    if( nl == 0 && o ==  1 ) return 2;
    if( nl == 0 && o == -1 ) return 1;
    if( nl == 1 && o ==  1 ) return 0;
    if( nl == 1 && o == -1 ) return 2;
    if( nl == 2 && o ==  1 ) return 1;
    if( nl == 2 && o == -1 ) return 0;
    assert(false); 
}

int MeshManifold2D::get_next_neighbor( int nl, int o )
{
    assert( 0 <= nl && nl < 3 );
    assert( o == 1 || o == -1 );
    if( nl == 0 && o ==  1 ) return 1;
    if( nl == 0 && o == -1 ) return 2;
    if( nl == 1 && o ==  1 ) return 2;
    if( nl == 1 && o == -1 ) return 0;
    if( nl == 2 && o ==  1 ) return 0;
    if( nl == 2 && o == -1 ) return 1;
    assert(false); 
}

int MeshManifold2D::get_prev_vertex  ( int vl, int o )
{
    assert( 0 <= vl && vl < 3 );
    assert( o == 1 || o == -1 );
    if( vl == 0 && o ==  1 ) return 2;
    if( vl == 0 && o == -1 ) return 1;
    if( vl == 1 && o ==  1 ) return 0;
    if( vl == 1 && o == -1 ) return 2;
    if( vl == 2 && o ==  1 ) return 1;
    if( vl == 2 && o == -1 ) return 0;
    assert(false); 
}

int MeshManifold2D::get_next_vertex  ( int vl, int o )
{
    assert( 0 <= vl && vl < 3 );
    assert( o == 1 || o == -1 );
    if( vl == 0 && o ==  1 ) return 1;
    if( vl == 0 && o == -1 ) return 2;
    if( vl == 1 && o ==  1 ) return 2;
    if( vl == 1 && o == -1 ) return 0;
    if( vl == 2 && o ==  1 ) return 0;
    if( vl == 2 && o == -1 ) return 1;
    assert(false); 
}

int MeshManifold2D::get_first_vertex( int o )
{
    assert( o == 1 || o == -1 );
    if( o == 1 ) return 0; else return 1; 
}

int MeshManifold2D::get_last_vertex ( int o )
{
    assert( o == 1 || o == -1 );
    if( o == 1 ) return 1; else return 0;
}


      








/* parents of an edge */

bool MeshManifold2D::is_triangle_edge_parent( int t, int e ) const
{
    assert( 0 <=  t &&  t < counter_triangles );
    assert( 0 <=  e &&  e < counter_edges );
    
    return ( data_edge_parents[e][0] == t ) || ( data_edge_parents[e][1] == t );
}


int MeshManifold2D::indexof_triangle_edge_parent( int t, int e ) const
{
    assert( 0 <=  t &&  t < counter_triangles );
    assert( 0 <=  e &&  e < counter_edges );
    
    if     ( data_edge_parents[e][0] == t ) return 0;
    else if( data_edge_parents[e][1] == t ) return 1;
    else                                    assert(false);
}


int MeshManifold2D::count_triangle_edge_parents( int e ) const 
{
    assert( 0 <=  e &&  e < counter_edges );
    
    if     ( data_edge_parents[e][0] != nullindex && data_edge_parents[e][1] != nullindex ) return 2;
    else if( data_edge_parents[e][0] != nullindex && data_edge_parents[e][1] == nullindex ) return 1;
    else if( data_edge_parents[e][0] == nullindex && data_edge_parents[e][1] != nullindex ) return 1;
    else                                                                                    assert(false);
}

const std::array<int,2> MeshManifold2D::get_triangle_edge_parents( int e ) const
{
    assert( 0 <=  e &&  e < counter_edges );
    
    return data_edge_parents[e];
}

int MeshManifold2D::get_edge_othertriangleparent( int t, int e ) const
{
    assert( 0 <=  t &&  t < counter_triangles );
    assert( 0 <=  e &&  e < counter_edges );
    
    const auto& parents = get_triangle_edge_parents(e);
    if( parents[0] == t ) return parents[1]; else return parents[0];
}

int MeshManifold2D::orientation_induced( int t, int el ) const
{
    assert( 0 <= t  && t  < counter_triangles );
    assert( 0 <= el && el < 3 );
    
    const std::array<int,2>& edge_global = get_edge_vertices( data_triangle_edges[t][el] );
    const std::array<int,2>& edge_local  = duple_from_triple( data_triangle_vertices[t], edgeindex_to_vertexindices( el ) );
    
    assert( vertexpairs_equivalent( edge_global, edge_local ) );
    
    if( edge_global[0] == edge_local[0] ) return 1; else return -1;
}




/* triangle parents of a vertex */

int MeshManifold2D::indexof_triangle_vertex_parent( int t, int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  std::vector<int> triangles = get_triangle_parents_of_vertex( v );
  
  auto iter = std::find( triangles.begin(), triangles.end(), t ); 
  assert( iter != triangles.end() );
  
  return iter - triangles.begin();
}

int MeshManifold2D::count_vertex_triangle_parents( int v ) const
{
  return get_triangle_parents_of_vertex( v ).size();
}

const std::vector<int> MeshManifold2D::get_triangle_parents_of_vertex( int v ) const
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

int MeshManifold2D::indexof_edge_vertex_parent( int e, int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  std::vector<int> edges = get_edge_parents_of_vertex( v );
  
  auto iter = std::find( edges.begin(), edges.end(), e );
  assert( iter != edges.end() );
  
  return iter - edges.begin();
}

int MeshManifold2D::count_vertex_edge_parents( int v ) const
{
  return get_edge_parents_of_vertex( v ).size();
}

const std::vector<int> MeshManifold2D::get_edge_parents_of_vertex( int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  
  std::vector<int> ret;
  
  for( int e = 0; e < count_edges(); e++ ) 
    for( int ev : get_edge_vertices(e) )
      if( v == ev )
        ret.push_back( e );
  
  return ret;
}




















/* index conversion */

// edge order : [01], [02], [12]
// the neigborindex and the edgeindex correspond to each other 
// opposing vertices: 2, 1, 0



int MeshManifold2D::neighborindex_to_edgeindex( int n )
{
    return n;
}

int MeshManifold2D::edgeindex_to_neighborindex( int e )
{
    return e;
}

std::array<int,2> MeshManifold2D::edgeindex_to_vertexindices( int e )
{
    if     ( e == 0 ) return {0,1};
    else if( e == 1 ) return {0,2};
    else if( e == 2 ) return {1,2};
    else assert(false);
}

std::array<int,2> MeshManifold2D::duple_from_triple( std::array<int,3> triple, std::array<int,2> duple )
{
    return { triple[duple[0]], triple[duple[1]] };
}

int MeshManifold2D::vertexindex_to_opposing_edgeindex( int v )
{
    if     ( v == 0 ) return 2;
    else if( v == 1 ) return 1;
    else if( v == 2 ) return 0;
    else assert(false);
}

int MeshManifold2D::edgeindex_to_opposing_vertexindex( int e )
{
    if     ( e == 0 ) return 2;
    else if( e == 1 ) return 1;
    else if( e == 2 ) return 0;
    else assert(false);
}

bool MeshManifold2D::vertexpairs_equivalent( std::array<int,2> e1, std::array<int,2> e2 )
{
    assert( e1[0] != e1[1] && e2[0] != e2[1] );
    
    if     ( e1[0] == e2[0] && e1[1] == e2[1] ) return true;
    else if( e1[0] == e2[1] && e1[1] == e2[0] ) return true;
    else                                        return false;
}








// // int MeshManifold2D::is_vertexparents_cyclic( int v );
// // int MeshManifold2D::count_vertexparents( int v, int* ts );
// // int MeshManifold2D::list_vertexparents( int v, int* ts );
// // int MeshManifold2D::count_vertexparents_edges( int v, int* es );
// // int MeshManifold2D::list_vertexparents_edges( int v, int* es );



/* bisect a single edge */

void MeshManifold2D::bisect_edge( int e )
{
    assert( 0 <= e && e < counter_edges );
    
    if( count_triangle_edge_parents( e ) == 1 )
        bisect_outer_edge( e );
    else if( count_triangle_edge_parents( e ) == 2 )
        bisect_inner_edge( e );
    else
        assert(false);
}

void MeshManifold2D::bisect_outer_edge( int e )
{
    std::cout << "Outer bisection" << nl;
  
    assert( 0 <= e && e < counter_edges );
    assert( count_triangle_edge_parents( e ) == 1 );
    
    int o = 1; /* TODO: Bisection of outer simplices abgleichen mit bisection innerer simplizes */
    
    FloatVector midpoint = get_edge_midpoint( e );
    getcoordinates().append( midpoint );
    
    int t            = get_triangle_edge_parents(e)[0];
    int e_local      = indexof_triangle_edge( t, e );
    int e_local_prev = get_prev_edge( e_local, o );
    int e_local_next = get_next_edge( e_local, o );
    
    assert( e_local != e_local_prev && e_local != e_local_next && e_local_next != e_local_prev );
    
    assert( data_triangle_edges[t][e_local] == get_triangle_edges(t)[e_local] );
    assert( e == get_triangle_edges(t)[e_local] );
    
    int e_prev = get_triangle_edges(t)[ e_local_prev ];
    int e_next = get_triangle_edges(t)[ e_local_next ];
    int n_prev = get_triangle_neighbors(t)[ edgeindex_to_neighborindex(e_local_prev) ];
    int n_next = get_triangle_neighbors(t)[ edgeindex_to_neighborindex(e_local_next) ];
    
    assert( e != e_prev && e != e_next && e_next != e_prev );
    assert( t != n_prev && t != n_next );
    if( n_prev != nullindex ) assert( n_next != n_prev );
    
    int v_local_opp  = edgeindex_to_opposing_vertexindex( e_local );
    int v_local_forw = edgeindex_to_opposing_vertexindex( e_local_prev );
    int v_local_back = edgeindex_to_opposing_vertexindex( e_local_next );
    
    assert( v_local_opp != v_local_forw && v_local_opp != v_local_back && v_local_forw != v_local_back );
    
    std::array<int,3> t_new_vertices_back 
      = { get_triangle_vertices(t)[v_local_opp ],
          get_triangle_vertices(t)[v_local_back],
          counter_vertices
        };
    
    std::array<int,3> t_new_vertices_forw
      = { counter_vertices,
          get_triangle_vertices(t)[v_local_forw],
          get_triangle_vertices(t)[v_local_opp ]
        };
    
    std::array<int,3> t_new_edges_back;
    t_new_edges_back[ vertexindex_to_opposing_edgeindex(0) ] = e;
    t_new_edges_back[ vertexindex_to_opposing_edgeindex(1) ] = counter_edges + 1;
    t_new_edges_back[ vertexindex_to_opposing_edgeindex(2) ] = e_prev;
    
    std::array<int,3> t_new_edges_forw;
    t_new_edges_forw[ vertexindex_to_opposing_edgeindex(0) ] = e_next;
    t_new_edges_forw[ vertexindex_to_opposing_edgeindex(1) ] = counter_edges + 1;
    t_new_edges_forw[ vertexindex_to_opposing_edgeindex(2) ] = counter_edges;
    
    std::array<int,2> e_new_parents_back = {                 t, nullindex };
    std::array<int,2> e_new_parents_forw = { counter_triangles, nullindex };
    std::array<int,2> e_new_pair_of_t    = { counter_triangles, t         };
    
    
    data_triangle_vertices[t] = t_new_vertices_back;
    data_triangle_vertices.push_back( t_new_vertices_forw );
    
    data_triangle_edges   [t] = t_new_edges_back;
    data_triangle_edges.push_back   ( t_new_edges_forw );
    
    data_edge_parents[e] = e_new_parents_back;
    data_edge_parents.push_back( e_new_parents_forw );
    data_edge_parents.push_back( e_new_pair_of_t );
    
//     std::cout << "test: " << e_local << space << e_local_prev << space << e_local_next << nl;
//     std::cout << "test: " << e << space << e_prev << space << e_next << nl;
//     std::cout << "test: " << t << space << n_prev << space << n_next << nl;
    data_edge_parents[ e_prev ] = {                 t, n_prev };
    data_edge_parents[ e_next ] = { counter_triangles, n_next };
    
    counter_triangles += 1;
    counter_edges     += 2;
    counter_vertices  += 1;
    
}

void MeshManifold2D::bisect_inner_edge( int e )
{
    std::cout << "Inner bisection" << nl;
  
    assert( 0 <= e && e < counter_edges );
    assert( count_triangle_edge_parents( e ) == 2 );
    
    FloatVector midpoint = get_edge_midpoint( e );
    getcoordinates().append( midpoint );
    
    /* gather auxiliary variables */
    
    /* TODO: extract vertices, use them to declare orientations */
    
    int t1 = get_triangle_edge_parents(e)[0];
    int t2 = get_triangle_edge_parents(e)[1];
    
    int edge_vertex_backward = get_edge_vertices( e )[0];
    int edge_vertex_forward  = get_edge_vertices( e )[1];
    
    int e1_local      = indexof_triangle_edge( t1, e );
    int e2_local      = indexof_triangle_edge( t2, e );
    
    int v1_local_opp  = edgeindex_to_opposing_vertexindex( e1_local );
    int v2_local_opp  = edgeindex_to_opposing_vertexindex( e2_local );
    
    int v1_local_back = indexof_triangle_vertex( t1, edge_vertex_backward );
    int v2_local_back = indexof_triangle_vertex( t2, edge_vertex_backward );
    int v1_local_forw = indexof_triangle_vertex( t1, edge_vertex_forward  );
    int v2_local_forw = indexof_triangle_vertex( t2, edge_vertex_forward  );
    
//     int v1_local_back = edgeindex_to_opposing_vertexindex( e1_local_prev );
//     int v2_local_forw = edgeindex_to_opposing_vertexindex( e2_local_prev );
//     int v1_local_back = edgeindex_to_opposing_vertexindex( e1_local_next );
//     int v2_local_back = edgeindex_to_opposing_vertexindex( e2_local_next );
    
    assert( v1_local_opp != v1_local_forw && v1_local_opp != v1_local_back && v1_local_forw != v1_local_back );
    assert( v2_local_opp != v2_local_forw && v2_local_opp != v2_local_back && v2_local_forw != v2_local_back );
    
    assert( data_triangle_vertices[t1][v1_local_back] == data_triangle_vertices[t2][v2_local_back] );
    assert( data_triangle_vertices[t1][v1_local_forw] == data_triangle_vertices[t2][v2_local_forw] );
    
//     int o1 = orientation_induced( t1, indexof_triangle_edge( t1, e ) );
//     int o2 = orientation_induced( t2, indexof_triangle_edge( t2, e ) );
    
    int e1_local_prev = vertexindex_to_opposing_edgeindex( v1_local_forw );
    int e2_local_prev = vertexindex_to_opposing_edgeindex( v2_local_forw );
    int e1_local_next = vertexindex_to_opposing_edgeindex( v1_local_back );
    int e2_local_next = vertexindex_to_opposing_edgeindex( v2_local_back );
    
//     int e1_local_prev = get_backwarding_edge( t1, e1_local, o1 );
//     int e2_local_prev = get_backwarding_edge( t2, e2_local, o2 );
//     int e1_local_next = get_forwarding_edge( t1, e1_local, o1 );
//     int e2_local_next = get_forwarding_edge( t2, e2_local, o2 );
    
    assert( e1_local != e1_local_prev && e1_local != e1_local_next && e1_local_next != e1_local_prev );
    assert( e2_local != e2_local_prev && e2_local != e2_local_next && e2_local_next != e2_local_prev );
    
    assert( data_triangle_edges[t1][e1_local] == get_triangle_edges(t1)[e1_local] );
    assert( data_triangle_edges[t2][e2_local] == get_triangle_edges(t2)[e2_local] );
    assert( e == get_triangle_edges(t1)[e1_local] );
    assert( e == get_triangle_edges(t2)[e2_local] );
    
    int e1_prev = get_triangle_edges(t1)[ e1_local_prev ];
    int e2_prev = get_triangle_edges(t2)[ e2_local_prev ];
    int e1_next = get_triangle_edges(t1)[ e1_local_next ];
    int e2_next = get_triangle_edges(t2)[ e2_local_next ];
    int n1_prev = get_triangle_neighbors(t1)[ edgeindex_to_neighborindex(e1_local_prev) ];
    int n2_prev = get_triangle_neighbors(t2)[ edgeindex_to_neighborindex(e2_local_prev) ];
    int n1_next = get_triangle_neighbors(t1)[ edgeindex_to_neighborindex(e1_local_next) ];
    int n2_next = get_triangle_neighbors(t2)[ edgeindex_to_neighborindex(e2_local_next) ];
    
    assert( e != e1_prev && e != e1_next && e1_next != e1_prev );
    assert( e != e2_prev && e != e2_next && e2_next != e2_prev );
    assert( t1 != n1_prev && t1 != n1_next );
    assert( t2 != n2_prev && t2 != n2_next );
    if( n1_prev != nullindex ) assert( n1_next != n1_prev );
    if( n2_prev != nullindex ) assert( n2_next != n2_prev );
    
    
    /* assemble new data */
    
    int t1_new = counter_triangles;
    int t2_new = counter_triangles + 1;
        
    int e_half_back = e;
    int e_half_forw = counter_edges;
    int e1_section  = counter_edges + 1;
    int e2_section  = counter_edges + 2;
    
    std::array<int,3> t1_new_vertices_back 
      = { get_triangle_vertices(t1)[v1_local_opp ],
          get_triangle_vertices(t1)[v1_local_back],
          counter_vertices
        };
    
    std::array<int,3> t2_new_vertices_back 
      = { get_triangle_vertices(t2)[v2_local_opp ],
          get_triangle_vertices(t2)[v2_local_back],
          counter_vertices
        };
    
    std::array<int,3> t1_new_vertices_forw
      = { counter_vertices,
          get_triangle_vertices(t1)[v1_local_forw],
          get_triangle_vertices(t1)[v1_local_opp ]
        };
    
    std::array<int,3> t2_new_vertices_forw
      = { counter_vertices,
          get_triangle_vertices(t2)[v2_local_forw],
          get_triangle_vertices(t2)[v2_local_opp ]
        };
        
    std::array<int,3> t1_new_edges_back;
    t1_new_edges_back[ vertexindex_to_opposing_edgeindex(0) ] = e_half_back;
    t1_new_edges_back[ vertexindex_to_opposing_edgeindex(1) ] = e1_section;
    t1_new_edges_back[ vertexindex_to_opposing_edgeindex(2) ] = e1_prev;
    
    std::array<int,3> t1_new_edges_forw;
    t1_new_edges_forw[ vertexindex_to_opposing_edgeindex(0) ] = e1_next;
    t1_new_edges_forw[ vertexindex_to_opposing_edgeindex(1) ] = e1_section;
    t1_new_edges_forw[ vertexindex_to_opposing_edgeindex(2) ] = e_half_forw;
    
    std::array<int,3> t2_new_edges_back;
    t2_new_edges_back[ vertexindex_to_opposing_edgeindex(0) ] = e_half_back;
    t2_new_edges_back[ vertexindex_to_opposing_edgeindex(1) ] = e2_section;
    t2_new_edges_back[ vertexindex_to_opposing_edgeindex(2) ] = e2_prev;
    
    std::array<int,3> t2_new_edges_forw;
    t2_new_edges_forw[ vertexindex_to_opposing_edgeindex(0) ] = e2_next;
    t2_new_edges_forw[ vertexindex_to_opposing_edgeindex(1) ] = e2_section;
    t2_new_edges_forw[ vertexindex_to_opposing_edgeindex(2) ] = e_half_forw;
    
    std::array<int,2> e_new_parents_back = { t1    , t2     };
    std::array<int,2> e_new_parents_forw = { t1_new, t2_new };
    std::array<int,2> e_new_pair_of_t1   = { t1    , t1_new };
    std::array<int,2> e_new_pair_of_t2   = { t2    , t2_new };
    
    /* transfer the new data */
    
    data_triangle_vertices[t1] = t1_new_vertices_back;
    data_triangle_vertices[t2] = t2_new_vertices_back;
    data_triangle_vertices.push_back( t1_new_vertices_forw );
    data_triangle_vertices.push_back( t2_new_vertices_forw );
    
    data_triangle_edges   [t1] = t1_new_edges_back;
    data_triangle_edges   [t2] = t2_new_edges_back;
    data_triangle_edges.push_back   ( t1_new_edges_forw );
    data_triangle_edges.push_back   ( t2_new_edges_forw );
    
    data_edge_parents[e] = e_new_parents_back;
    data_edge_parents.push_back( e_new_parents_forw );
    data_edge_parents.push_back( e_new_pair_of_t1 );
    data_edge_parents.push_back( e_new_pair_of_t2 );
    
    data_edge_parents[ e1_prev ] = { t1    , n1_prev };
    data_edge_parents[ e2_prev ] = { t2    , n2_prev };
    data_edge_parents[ e1_next ] = { t1_new, n1_next };
    data_edge_parents[ e2_next ] = { t2_new, n2_next };
    
    counter_triangles += 2;
    counter_edges     += 3;
    counter_vertices  += 1;
}





void MeshManifold2D::uniformrefinement()
{
    int offset_triangles = counter_triangles;
    int offset_edges     = counter_edges;
    int offset_vertices  = counter_vertices;
    
    /* Add new vertices */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      const auto mid = get_edge_midpoint( e );
      getcoordinates().append( mid );
    }
    
    
    /* update the list of triangles -> vertices */
    
    data_triangle_vertices.resize( 4 * counter_triangles );
    
    
    /* sweep over the triangles and update/create the data */
    
    for( int t = 0; t < counter_triangles; t++ )
    {
        /* Gather auxiliary data */
      
        int v_00 = data_triangle_vertices[t][0];
        int v_11 = data_triangle_vertices[t][1];
        int v_22 = data_triangle_vertices[t][2];
        
        int v_01 = offset_vertices + data_triangle_edges[t][ 0 ];
        int v_02 = offset_vertices + data_triangle_edges[t][ 1 ];
        int v_12 = offset_vertices + data_triangle_edges[t][ 2 ];
        
        std::array<int,3> new_triangle_vertices_0, new_triangle_vertices_1, new_triangle_vertices_2, new_triangle_vertices_m;
        
        new_triangle_vertices_0[0] = v_00;
        new_triangle_vertices_0[1] = v_01;
        new_triangle_vertices_0[2] = v_02;
        
        new_triangle_vertices_1[0] = v_11;
        new_triangle_vertices_1[1] = v_01;
        new_triangle_vertices_1[2] = v_12;
        
        new_triangle_vertices_2[0] = v_22;
        new_triangle_vertices_2[1] = v_02;
        new_triangle_vertices_2[2] = v_12;
        
        new_triangle_vertices_m[0] = v_01;
        new_triangle_vertices_m[1] = v_02;
        new_triangle_vertices_m[2] = v_12;
        
        data_triangle_vertices[ t + 0 * offset_triangles ] = new_triangle_vertices_0;
        data_triangle_vertices[ t + 1 * offset_triangles ] = new_triangle_vertices_1;
        data_triangle_vertices[ t + 2 * offset_triangles ] = new_triangle_vertices_2;
        data_triangle_vertices[ t + 3 * offset_triangles ] = new_triangle_vertices_m;
        
    }
    
    
    
    /* Count edges, T->E, E->T  */
    
    /* delete all edges and fill up with nullindices */ 
    
    data_triangle_edges.resize(0);
    data_triangle_edges.resize( 4 * counter_triangles, { nullindex, nullindex, nullindex } );
    
    /* delete all edge parent index */ 
    
    data_edge_parents.resize(0);
    data_edge_parents.reserve( 2 * counter_edges + 3 * counter_triangles );
    
    /* first the inner edges */
    
    for( int t1 = 0; t1 < data_triangle_vertices.size(); t1++ )
    for( int t2 = 0; t2 < t1                           ; t2++ )
    for( int e1 = 0; e1 < 3; e1++ )
    for( int e2 = 0; e2 < 3; e2++ )
    {
        const auto& ve1 = duple_from_triple( data_triangle_vertices[t1], edgeindex_to_vertexindices( e1 ) );
        const auto& ve2 = duple_from_triple( data_triangle_vertices[t2], edgeindex_to_vertexindices( e2 ) );
        if( vertexpairs_equivalent( ve1, ve2 ) )
        {
            data_edge_parents.push_back( { t1, t2 } );
            data_triangle_edges[t1][ e1 ] = data_edge_parents.size() - 1;
            data_triangle_edges[t2][ e2 ] = data_edge_parents.size() - 1;
        }
    }
    
    /* then the outer edges */
    
    for( int t = 0; t < data_triangle_vertices.size(); t++ )
    for( int e = 0; e < 3; e++ )
    {
        if( data_triangle_edges[t][e] == nullindex ) 
        {
            data_edge_parents.push_back( { t, nullindex } );
            data_triangle_edges[t][e] = data_edge_parents.size() - 1;
        }
    }
    
//     counter_edges = data_edge_parents.size();
    
    /* update the counters */
    
    counter_vertices  = counter_vertices + counter_edges;
    counter_edges     = 2 * counter_edges + 3 * counter_triangles;
    counter_triangles = 4 * counter_triangles;
    
}




// void MeshManifold2D::uniformrefinement()
// {
//     assert( false );
//     
//     int offset_t = counter_triangles;
//     int offset_e = counter_edges;
//     int offset_v = counter_vertices;
//     
//     /* Add new vertices */
//     
//     for( int e = 0; e < counter_edges; e++ )
//     {
//       const auto mid = get_edge_midpoint( e );
//       getcoordinates().append( mid );
//     }
//     
//     /* resize container */
//     
//     data_edges_parents.resize( 2 * counter_edges + 3 * counter_triangles );
//     
//     data_triangle_edges.resize( 4 * counter_triangles );
//     
//     data_triangle_vertices.resize( 4 * counter_triangles );
//     
//     
//     /* sweep over the triangles and update/create the data */
//     
//     /* TODO: We need to keep the old data throughout the whole process
//     /*       in order to keep up with the parents. Maybe there are other alternatives? */
//     
//     for( int t = 0; t < counter_triangles; t++ )
//     {
//         /* Gather auxiliary data */
//       
//         /* TODO: This is wrong ... */
//       
//         int e_01_v0 = data_triangle_edges[t][0];
//         int e_02_v0 = data_triangle_edges[t][1];
//         int e_12_v1 = data_triangle_edges[t][2];
//         
//         int e_01_v1 = counter_edges + data_triangle_edges[t][0];
//         int e_02_v2 = counter_edges + data_triangle_edges[t][1];
//         int e_12_v2 = counter_edges + data_triangle_edges[t][2];
//         
//         /* TODO: ... until here */
//       
//         int e_m0    = 2 * counter_edges + 0 * counter_triangles + t;
//         int e_m1    = 2 * counter_edges + 1 * counter_triangles + t;
//         int e_m2    = 2 * counter_edges + 2 * counter_triangles + t;
//         
//         int v_00 = data_triangle_vertices[t][0];
//         int v_11 = data_triangle_vertices[t][1];
//         int v_22 = data_triangle_vertices[t][2];
//         
//         int v_01 = offset_vertices + offset_edges * vertexindex_to_opposing_edgeindex(2)
//                    + data_triangle_edges[t][ vertexindex_to_opposing_edgeindex(2) ];
//         int v_02 = offset_vertices + offset_edges * vertexindex_to_opposing_edgeindex(1) 
//                    + data_triangle_edges[t][ vertexindex_to_opposing_edgeindex(1) ];
//         int v_12 = offset_vertices + offset_edges * vertexindex_to_opposing_edgeindex(0) 
//                    + data_triangle_edges[t][ vertexindex_to_opposing_edgeindex(0) ];
//         
//         /* TODO: gather auxiliary data about parents */
//         
//         /* Prepare new data sets */
//       
//         std:array<int,3> new_triangle_edges_0,    new_triangle_edges_1,    new_triangle_edges_2,    new_triangle_edges_m;
//         std:array<int,3> new_triangle_vertices_0, new_triangle_vertices_1, new_triangle_vertices_2, new_triangle_vertices_m;
//         
//         new_triangle_vertices_0[0] = v_00;      new_triangle_edges_0[0] = ;
//         new_triangle_vertices_0[1] = v_01;      new_triangle_edges_0[1] = ;
//         new_triangle_vertices_0[2] = v_02;      new_triangle_edges_0[2] = ;
//         
//         new_triangle_vertices_1[0] = v_11;      new_triangle_edges_1[0] = ;
//         new_triangle_vertices_1[1] = v_01;      new_triangle_edges_1[1] = ;
//         new_triangle_vertices_1[2] = v_12;      new_triangle_edges_1[2] = ;
//         
//         new_triangle_vertices_2[0] = v_22;      new_triangle_edges_2[0] = ;
//         new_triangle_vertices_2[1] = v_02;      new_triangle_edges_2[1] = ;
//         new_triangle_vertices_2[2] = v_12;      new_triangle_edges_2[2] = ;
//         
//         new_triangle_vertices_m[0] = v_01;      new_triangle_edges_m[0] = ;
//         new_triangle_vertices_m[1] = v_02;      new_triangle_edges_m[1] = ;
//         new_triangle_vertices_m[2] = v_12;      new_triangle_edges_m[2] = ;
//         
//         /* Update data */
//       
//         
//     }
//     
//     
//     /* update the counters */
//     
//     counter_vertices  = counter_edges;
//     counter_edges     = 2 * counter_edges + 3 * counter_triangles;
//     counter_triangles = 4 * counter_triangles
//     
// }
        








FloatVector MeshManifold2D::get_triangle_midpoint( int t )
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

FloatVector MeshManifold2D::get_edge_midpoint    ( int e )
{
    assert( 0 <= e && e < counter_edges );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = ( getcoordinates().getdata( get_edge_vertices(e)[0], d ) 
                 + getcoordinates().getdata( get_edge_vertices(e)[1], d )
               ) / 2.;
    return mid;
}


        



