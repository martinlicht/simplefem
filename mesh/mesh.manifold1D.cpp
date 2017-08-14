
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

#include "mesh.manifold1D.hpp"




MeshManifold1D::MeshManifold1D( int outerdim )
:
    Mesh( 1, outerdim ),
    
    counter_edges(0),
    counter_vertices(0),
    
    data_edge_vertices(0),
    data_vertex_firstparent(0),
    data_edge_nextparents(0)
{
    check();
}


MeshManifold1D::MeshManifold1D( 
    int outerdim,
    const Coordinates& coords,
    const std::vector<std::array<int,2>> edge_vertices
)
:
    Mesh( 1, outerdim ),
    
    counter_edges( edge_vertices.size() ),
    counter_vertices(0),
    
    data_edge_vertices( edge_vertices ),
    data_vertex_firstparent( 0 ),
    data_edge_nextparents( counter_edges, { nullindex, nullindex } )
{
    
    getcoordinates() = coords;
    
    /* 1. Count edges, transfer data */ 
    /* DONE */
    
    
    /* 2. Count vertices, allocate memory */
    counter_vertices = 0;
    for( const auto& duple : data_edge_vertices )
    for( const int& vertex : duple )
      counter_vertices = counter_vertices < vertex ? vertex : counter_vertices; 
    counter_vertices += 1;
    
    data_vertex_firstparent.resize( counter_vertices, nullindex );
    
    /* 3. For each vertex, set the first parent and the neighboring parents */
    
    for( int e =  0; e  <  counter_edges; e++  )
    for( int vi = 0; vi <=             1; vi++ )
    {
      int v = data_edge_vertices[e][vi];
      
      assert( 0 <= v && v < counter_vertices );
      
      if( data_vertex_firstparent[v] == nullindex ) {
        
        data_vertex_firstparent[v] = e;
        
      } else {
        
        int old_first_parent = data_vertex_firstparent[v];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_edges );
        assert( data_edge_nextparents[ e ][ vi ] == nullindex );
        
        data_vertex_firstparent[v] = e;
        data_edge_nextparents[ e ][ vi ] = old_first_parent;
        
      }
      
      assert( data_vertex_firstparent[v] != nullindex );
      assert( 0 <= data_vertex_firstparent[v] && data_vertex_firstparent[v] < counter_edges );
      
      
      
    }
    
    check();
}


MeshManifold1D::MeshManifold1D( 
    int outerdim,
    const Coordinates& coords,
    const std::vector<std::array<int,2>> edge_vertices,
    const std::vector<int              > vertex_firstparent,
    const std::vector<std::array<int,2>> edge_nextparents
)
:
    Mesh( 1, outerdim ),
    
    counter_edges( edge_vertices.size() ),
    counter_vertices( vertex_firstparent.size() ),
    
    data_edge_vertices( edge_vertices ),
    data_vertex_firstparent( vertex_firstparent ),
    data_edge_nextparents( edge_nextparents )
{
    
    getcoordinates() = coords;
    
    check();
}


MeshManifold1D::~MeshManifold1D()
{
    
}


void MeshManifold1D::check() const
{
    
    /* 1. Check the array sizes */
    assert( counter_edges == data_edge_vertices.size() );
    assert( counter_edges == data_edge_nextparents.size() );
    assert( counter_vertices == data_vertex_firstparent.size() );
    assert( count_vertices() == getcoordinates().getnumber() );
    
    
    /* * * * * Data integrity
     * 
     * each edge: each vertex is a valid index
     * each edge: each vertex is unique 
     * each edge: the next parents make sense 
     * 
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
        assert( data_edge_vertices[e][0] != nullindex );
        assert( data_edge_vertices[e][1] != nullindex );
        assert( data_edge_vertices[e][0] != data_edge_vertices[e][1] );
        
        if( data_edge_nextparents[e][0] != nullindex || data_edge_nextparents[e][1] != nullindex )
          assert( data_edge_nextparents[e][0] != data_edge_nextparents[e][1] );
        
        if( data_edge_nextparents[e][0] != nullindex )
          assert( 0 <= data_edge_nextparents[e][0] && data_edge_nextparents[e][0] < counter_edges );
        
        if( data_edge_nextparents[e][0] != nullindex )
          assert( data_edge_vertices[ data_edge_nextparents[e][0] ][0] == data_edge_vertices[e][0] 
                  ||
                  data_edge_vertices[ data_edge_nextparents[e][0] ][1] == data_edge_vertices[e][0] );
        
        if( data_edge_nextparents[e][1] != nullindex )
          assert( 0 <= data_edge_nextparents[e][1] && data_edge_nextparents[e][1] < counter_edges );

        if( data_edge_nextparents[e][1] != nullindex )
          assert( data_edge_vertices[ data_edge_nextparents[e][1] ][0] == data_edge_vertices[e][1] 
                  ||
                  data_edge_vertices[ data_edge_nextparents[e][1] ][1] == data_edge_vertices[e][1] );
    }
    
    /* * * * * * Data integrity
     * 
     * each vertex_firstparent: first parent is non-null 
     * 
     */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
        int p = data_vertex_firstparent[v];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_edges );
        
        assert( data_edge_vertices[p][0] == v || data_edge_vertices[p][1] == v );
    }
    
    /* * * * * * Data integrity
     * 
     * check that each is listed as parent somewhere 
     * 
     */
    
    for( int e  = 0; e  <  counter_edges; e++ )
    for( int vi = 0; vi <=             1; vi++ )
    {
      
      int v = data_edge_vertices[e][vi];
      
      int p = data_vertex_firstparent[v];
      
      assert( p != nullindex );
      
      while( p != e && p != nullindex )
        p = data_edge_nextparents[p][ ( data_edge_vertices[p][0] == v ) ? 0 : 1 ];
      
      assert( p == e );
      
    }
    
  
}






void MeshManifold1D::print( std::ostream& os ) const
{
    os << "Printe Triangulation of 1D Manifold!" << std::endl;
    
    os << counter_edges << space << counter_vertices << nl;
    
    os << "Edge vertices" << std::endl;
    
    for( const auto& duple : data_edge_vertices )
      std::cout << duple[0] << space << duple[1] << nl;
    
    os << "Vertex first parents" << std::endl;
    
    for( int fp : data_vertex_firstparent )
      std::cout << fp << nl;
    
    os << "Edge next parents " << std::endl;
    
    for( const auto& duple : data_edge_nextparents )
      std::cout << duple[0] << space << duple[1] << nl;
    
    os << "Finished printing" << nl;
    
}






bool MeshManifold1D::dimensioncounted( int dim ) const
{
    assert( 0 <= dim && dim <= 1 );
    return true;
}

int MeshManifold1D::countsimplices( int dim ) const
{
  if( dim == 0 )
    return count_vertices();
  else if( dim == 1 )
    return count_edges();
  else
    assert(false);
}

bool MeshManifold1D::subsimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub < sup && sup <= 1 );
    return true;
}

const IndexMap MeshManifold1D::getsubsimplices( int sup, int sub, int cell ) const
{
  assert( 1 == sup );
  assert( 0 == sub );
  assert( 0 <= cell && cell < count_edges() );
  
  auto temp = get_edge_vertices(cell);
  return IndexMap( IndexRange(0,1), std::vector<int>( temp.begin(), temp.end() ) );
    
}

bool MeshManifold1D::supersimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub < sup && sup <= 1 );
    return true;
}

const std::vector<int> MeshManifold1D::getsupersimplices( int sup, int sub, int cell ) const
{
  assert( 1 == sup );
  assert( 0 == sub );
  assert( 0 <= cell && cell < count_vertices() );
  
  auto temp = get_edge_parents_of_vertex( cell ); 
  return std::vector<int>( temp.begin(), temp.end() );
}








/* Count number of elements */

int MeshManifold1D::count_edges() const
{
    return counter_edges;
}

int MeshManifold1D::count_vertices() const
{
    return counter_vertices;
}



/* subsimplex relation of edges and vertices */

bool MeshManifold1D::contains_edge_vertex( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    
    return ( data_edge_vertices[e][0] == v ) || ( data_edge_vertices[e][1] == v );
} 

int MeshManifold1D::indexof_edge_vertex( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    if     ( data_edge_vertices[e][0] == v ) return 0;
    else if( data_edge_vertices[e][1] == v ) return 1;
    else                                     assert(false);
} 

const std::array<int,2> MeshManifold1D::get_edge_vertices( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    return data_edge_vertices[e];
} 




/* edge parents of a vertex */

int MeshManifold1D::count_vertex_edge_parents( int v ) const
{
  return get_edge_parents_of_vertex( v ).size();
}

int MeshManifold1D::get_vertex_firstparent( int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  return data_vertex_firstparent[ v ];
}

int MeshManifold1D::get_vertex_nextparent( int v, int e ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= e && e < counter_edges    );
  
  if( data_edge_vertices[e][0] == v )
    return data_edge_nextparents[e][0];
  else if( data_edge_vertices[e][1] == v )
    return data_edge_nextparents[e][1];
  else
    assert(false);
}

bool MeshManifold1D::is_edge_vertex_parent( int e, int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= e && e < counter_edges    );
  return data_edge_vertices[e][0] == v || data_edge_vertices[e][1] == v;
}

int MeshManifold1D::indexof_edge_vertex_parent( int e, int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  std::vector<int> edges = get_edge_parents_of_vertex( v );
  
  auto iter = std::find( edges.begin(), edges.end(), e ); 
  assert( iter != edges.end() );
  
  return iter - edges.begin();
}

std::vector<int> MeshManifold1D::get_edge_parents_of_vertex( int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  
  std::vector<int> ret;
  
  for( int e = 0; e < count_edges(); e++ ) 
    for( int ev : get_edge_vertices(e) )
      if( v == ev )
        ret.push_back( e );
  
  return ret;
}



void MeshManifold1D::bisect_edge( int e )
{
    assert( 0 <= e && e < counter_edges );
    check();
    
    /* Collect the old data */
    
    int vertex_back  = data_edge_vertices[e][0];
    int vertex_front = data_edge_vertices[e][1];
    int nextparent_back  = data_edge_nextparents[e][0];
    int nextparent_front = data_edge_nextparents[e][1];
    int firstparent_back  = data_vertex_firstparent[vertex_back ];
    int firstparent_front = data_vertex_firstparent[vertex_front];
    
    int back_previousparent   = nullindex;
    int back_previousparent_localindex = nullindex;
    if( e != get_vertex_firstparent( vertex_back ) )
      for( back_previousparent = get_vertex_firstparent( vertex_back  );
           back_previousparent != nullindex && get_vertex_nextparent(vertex_back,back_previousparent) != e; 
           back_previousparent = get_vertex_nextparent(vertex_back,back_previousparent) 
         ); 
    if( back_previousparent  != nullindex ) 
      back_previousparent_localindex  = indexof_edge_vertex( back_previousparent,  vertex_back  );
    
    int front_previousparent  = nullindex;
    int front_previousparent_localindex = nullindex;
    if( e != get_vertex_firstparent( vertex_front ) )
      for( front_previousparent = get_vertex_firstparent( vertex_front );
           front_previousparent != nullindex && get_vertex_nextparent(vertex_front,front_previousparent) != e; 
           front_previousparent = get_vertex_nextparent(vertex_front,front_previousparent) 
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
    
    int back_backnextparent   = nextparent_back;
    int back_frontnextparent  = ne;
    int front_backnextparent  = nullindex;
    int front_frontnextparent = nextparent_front;
    
    int firstparent_newvertex = nv;
    
    /* Allocate memory */
    
    data_edge_nextparents.resize  ( counter_edges    + 1 );
    data_edge_vertices.resize     ( counter_edges    + 1 );
    data_vertex_firstparent.resize( counter_vertices + 1 );
    
    /* Write in the data */
    
    data_edge_vertices[e ][0] = back_backvertex;
    data_edge_vertices[e ][1] = back_frontvertex;
    data_edge_vertices[ne][0] = front_backvertex;
    data_edge_vertices[ne][1] = front_frontvertex;
    
    data_edge_nextparents[e ][0] = back_backnextparent;
    data_edge_nextparents[e ][1] = back_frontnextparent;
    data_edge_nextparents[ne][0] = front_backnextparent;
    data_edge_nextparents[ne][1] = front_frontnextparent;
    
    data_vertex_firstparent[nv] = e;
    
    if( back_previousparent  != nullindex ) {
      assert( data_edge_nextparents[ back_previousparent ][ back_previousparent_localindex ] == e );
      data_edge_nextparents[ back_previousparent ][ back_previousparent_localindex ] = e;
    } else {
      assert( data_vertex_firstparent[ vertex_back ] == e );
      data_vertex_firstparent[ vertex_back ] = e;
    }
    
    if( front_previousparent != nullindex ) {
      assert( data_edge_nextparents[ front_previousparent ][ front_previousparent_localindex ] == ne );
      data_edge_nextparents[ front_previousparent ][ front_previousparent_localindex ] = ne;
    } else {
      assert( data_vertex_firstparent[ vertex_front ] == e );
      data_vertex_firstparent[ vertex_front ] = ne;
    }    
    
    getcoordinates().append( midcoordinate );
        
    /* Update counter */
    counter_edges++;
    counter_vertices++;
    
    /* Done */
    
    
    check();
    
}


void MeshManifold1D::uniformrefinement()
{
    int old_counter_edges    = counter_edges;
    int old_counter_vertices = counter_vertices;
    
    check();
    
    data_edge_nextparents.reserve  ( 2 * old_counter_edges );
    data_edge_vertices.reserve     ( 2 * old_counter_edges );
    data_vertex_firstparent.reserve( old_counter_vertices  );
    getcoordinates().addcapacity   ( old_counter_edges     );
    
    for( int e = 0; e < old_counter_edges; e++ )
      bisect_edge( e );
      
    check();
}



void MeshManifold1D::improved_uniformrefinement()
{
    check();
    
    /* resize the arrays */
    
    data_edge_nextparents.resize  ( counter_edges * 2 );
    data_edge_vertices.resize     ( counter_edges * 2 );
    data_vertex_firstparent.resize( counter_edges + counter_vertices );
    getcoordinates().addcoordinates( counter_edges );
    
    
    /* create the new coordinates and fill them up */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      getcoordinates().loadvector( e , get_edge_midpoint( e ) );
    }
    
    
    /* for each vertex, set the new parent edge */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      int p = data_vertex_firstparent[v];
      
      assert( p != nullindex );
      
      int vi = ( ( data_edge_vertices[p][0] == v ) ? 0 : 1 );
      
      assert( data_edge_vertices[p][0] == v || data_edge_vertices[p][1] == v );
      assert( data_edge_vertices[p][vi] == v );
      
      data_vertex_firstparent[v] = p + vi * counter_edges;
    }
    
    
    /* for each edge, relocate the data of the old vertices' parents */
    
    for( int e  = 0; e  < counter_edges;  e++ )
    for( int vi = 0; vi <=            1; vi++ )
    {
      int q = data_edge_nextparents[e][vi];
      
      int v = data_edge_vertices[e][vi];
      
      if( q == nullindex ) {
        
        data_edge_nextparents[e + vi * counter_edges ][vi] = nullindex;
        
      } else if( q != nullindex ) {
        
        int vinp = ( data_edge_vertices[q][0] == v ? 0 : 1 );
        
        assert( data_edge_vertices[q][0] == v || data_edge_vertices[q][1] == v );
        assert( data_edge_vertices[q][vinp] == v );
        
        data_edge_nextparents[ e + vi * counter_edges ][vi] = q + vinp * counter_edges;
      
      } 
      
    }
    
    
    
    /* for each new vertex, create the first and second parent */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      data_vertex_firstparent[counter_vertices + e] = e;
      
      data_edge_nextparents[e                ][1] = e + counter_edges;
      data_edge_nextparents[e + counter_edges][0] = nullindex;
    }
    
    
    /* for each edge, create the new vertices */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int vertex_back  = data_edge_vertices[e][0];
      int vertex_front = data_edge_vertices[e][1];
      
      data_edge_vertices[e                ][0] = vertex_back;
      data_edge_vertices[e                ][1] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][0] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][1] = vertex_front;
    }
    
    
    /* update the counters */
    
    counter_vertices += counter_edges;
    counter_edges += counter_edges;
    
    
    /* DONE */
    
    check();
}






FloatVector MeshManifold1D::get_edge_midpoint    ( int e )
{
    assert( 0 <= e && e < counter_edges );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = (   getcoordinates().getdata( get_edge_vertices(e)[0], d ) 
                 + getcoordinates().getdata( get_edge_vertices(e)[1], d )
                ) / 2.;
    return mid;
}


        



