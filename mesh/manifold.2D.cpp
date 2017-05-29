
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
    coordinates(coords),
    
    counter_triangles( triangle_vertex_list.size() ),
    counter_edges(0),
    counter_vertices(0),
    
    data_triangle_vertices( triangle_vertex_list ),
    data_triangle_edges( triangle_vertex_list.size(), {nullindex,nullindex,nullindex} ),
    data_edge_parents(0)
//     data_vertex_firstparent(0)
{
    
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
        if( vertexlists_equivalent( ve1, ve2 ) )
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


ManifoldTriangulation2D::~ManifoldTriangulation2D()
{
    
}


void ManifoldTriangulation2D::check() const
{
    /* * Data integrity
     * 
     * - each triangle: each vertex is non-null
     * - each triangle: each vertex is unique
     * - each triangle: each edge is non-null
     * - each triangle: each edge is unique
     * 
     */
    
    assert( counter_triangles == data_triangle_vertices.size() );
    assert( counter_triangles == data_triangle_edges.size() );
    
//     cout << count_vertices() << space << coordinates.getnumber() << nl;
//     assert( count_vertices() == coordinates.getnumber() );
    
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
        
//         std::cout << "vergleich: " << data_edge_parents[e][0] << space << data_edge_parents[e][1] << nl;
        
        assert( data_edge_parents[e][0] != data_edge_parents[e][1] );
        assert( data_edge_parents[e][0] != nullindex || data_edge_parents[e][1] != nullindex );
        assert( data_edge_parents[e][0] != nullindex );
        
        int p0 = data_edge_parents[e][0];
        
        assert( data_triangle_edges[p0][0] == e || data_triangle_edges[p0][1] == e || data_triangle_edges[p0][2] == e );
        
        int p1 = data_edge_parents[e][1];
        
        if( p1 == nullindex) continue; 
        
        assert( data_triangle_edges[p1][0] == e || data_triangle_edges[p1][1] == e || data_triangle_edges[p1][2] == e );
        
    }
    
    
    
}

void ManifoldTriangulation2D::print( std::ostream& os ) const
{
    os << "Printe Triangulation of 2D Manifold!" << std::endl;
    
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
      std::cout << get_edge_vertices(e)[0] << space << get_edge_vertices(e)[1] << nl;
    
    
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
    assert( counter_triangles == data_triangle_edges.size() );
    assert( counter_triangles == data_triangle_vertices.size() );
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
    assert( is_not_nullindex( t ) );
    assert( contains_triangle_edge( t, e ) );
    
    int ei = indexof_triangle_edge( t, e );
    std::array<int,2> vi = edgeindex_to_vertexindices( ei );
    
    return    ( data_triangle_vertices[t][vi[0]] == v )
           || ( data_triangle_vertices[t][vi[1]] == v )
           ;
} 







int ManifoldTriangulation2D::indexof_triangle_edge   ( int t, int e ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= e && e < counter_edges );
    if     ( data_triangle_edges[t][0] == e ) return 0;
    else if( data_triangle_edges[t][1] == e ) return 1;
    else if( data_triangle_edges[t][2] == e ) return 2;
    else                                      assert(false);
} 

int ManifoldTriangulation2D::indexof_triangle_vertex( int t, int v ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= v && v < counter_vertices );
    if     ( data_triangle_vertices[t][0] == v ) return 0;
    else if( data_triangle_vertices[t][1] == v ) return 1;
    else if( data_triangle_vertices[t][2] == v ) return 2;
    else                                         assert(false);
} 

int ManifoldTriangulation2D::indexof_edge_vertex    ( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    
    int t = get_edge_parents( e )[0];
    assert( contains_triangle_edge( t, e ) );
    
    int ei = indexof_triangle_edge( t, e );
    std::array<int,2> vi = edgeindex_to_vertexindices( ei );
    
    if     ( data_triangle_vertices[t][vi[0]] == v ) return vi[0];
    else if( data_triangle_vertices[t][vi[1]] == v ) return vi[1];
    else                                             assert(false); 
} 










const std::array<int,3> ManifoldTriangulation2D::get_triangle_edges   ( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    return data_triangle_edges[t];
} 

const std::array<int,3> ManifoldTriangulation2D::get_triangle_vertices( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    return data_triangle_vertices[t];
}

const std::array<int,2> ManifoldTriangulation2D::get_edge_vertices    ( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    
    int t = get_edge_parents( e )[0];
    assert( contains_triangle_edge( t, e ) );
    
    int el = indexof_triangle_edge( t, e );
    
    std::array<int,2> vi = edgeindex_to_vertexindices( el );
    
    return { data_triangle_vertices[t][ vi[0] ], data_triangle_vertices[t][ vi[1] ] };
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
        
        if( data_edge_parents[e][0] == t )
          ret[ei] = data_edge_parents[e][1];
        else if( data_edge_parents[e][1] == t )
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
    
    if( ! is_triangle_neighbor( t1, t2 )  ) return false;
    if( ! contains_triangle_edge( t1, e ) ) return false;
    if( ! contains_triangle_edge( t2, e ) ) return false;
    
    return true;
}

int ManifoldTriangulation2D::get_edge_between( int t1, int t2 ) const
{
    assert( 0 <= t1 && t1 < counter_triangles );
    assert( 0 <= t2 && t2 < counter_triangles );
    
    assert( is_triangle_neighbor( t1, t2 ) );
    
    int n = indexof_triangle_neighbor( t1, t2 );
    
    return neighborindex_to_edgeindex( n );
}








int ManifoldTriangulation2D::get_prev_edge    ( int el, int o )
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

int ManifoldTriangulation2D::get_next_edge    ( int el, int o )
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

int ManifoldTriangulation2D::get_prev_neighbor( int nl, int o )
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

int ManifoldTriangulation2D::get_next_neighbor( int nl, int o )
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

int ManifoldTriangulation2D::get_prev_vertex  ( int vl, int o )
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

int ManifoldTriangulation2D::get_next_vertex  ( int vl, int o )
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

int ManifoldTriangulation2D::get_first_vertex( int o )
{
    assert( o == 1 || o == -1 );
    if( o == 1 ) return 0; else return 1; 
}

int ManifoldTriangulation2D::get_last_vertex ( int o )
{
    assert( o == 1 || o == -1 );
    if( o == 1 ) return 1; else return 0;
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
    
    if     ( data_edge_parents[e][0] == t ) return 0;
    else if( data_edge_parents[e][1] == t ) return 1;
    else                                    assert(false);
}


int ManifoldTriangulation2D::count_edge_parents( int e ) const 
{
    assert( 0 <=  e &&  e < counter_edges );
    
    if     ( data_edge_parents[e][0] != nullindex && data_edge_parents[e][1] != nullindex ) return 2;
    else if( data_edge_parents[e][0] != nullindex && data_edge_parents[e][1] == nullindex ) return 1;
    else if( data_edge_parents[e][0] == nullindex && data_edge_parents[e][1] != nullindex ) return 1;
    else                                                                                    assert(false);
}

const std::array<int,2> ManifoldTriangulation2D::get_edge_parents( int e ) const
{
    assert( 0 <=  e &&  e < counter_edges );
    
    return data_edge_parents[e];
}

int ManifoldTriangulation2D::get_edge_otherparent( int t, int e ) const
{
    assert( 0 <=  t &&  t < counter_triangles );
    assert( 0 <=  e &&  e < counter_edges );
    
    const auto& parents = get_edge_parents(e);
    if( parents[0] == t ) return parents[1]; else return parents[0];
}

int ManifoldTriangulation2D::orientation_induced( int t, int el ) const
{
    assert( 0 <= t  && t  < counter_triangles );
    assert( 0 <= el && el < 3 );
    
    const std::array<int,2>& edge_global = get_edge_vertices( data_triangle_edges[t][el] );
    const std::array<int,2>& edge_local  = duple_from_triple( data_triangle_vertices[t], edgeindex_to_vertexindices( el ) );
    
    assert( vertexlists_equivalent( edge_global, edge_local ) );
    
    if( edge_global[0] == edge_local[0] ) return 1; else return -1;
}





/* index conversion */

// edge order : [01], [02], [12]
// the neigborindex and the edgeindex correspond to each other 
// opposing vertices: 2, 1, 0



int ManifoldTriangulation2D::neighborindex_to_edgeindex( int n )
{
    return n;
}

int ManifoldTriangulation2D::edgeindex_to_neighborindex( int e )
{
    return e;
}

std::array<int,2> ManifoldTriangulation2D::edgeindex_to_vertexindices( int e )
{
    if     ( e == 0 ) return {0,1};
    else if( e == 1 ) return {0,2};
    else if( e == 2 ) return {1,2};
    else assert(false);
}

std::array<int,2> ManifoldTriangulation2D::duple_from_triple( std::array<int,3> triple, std::array<int,2> duple )
{
    return { triple[duple[0]], triple[duple[1]] };
}

int ManifoldTriangulation2D::vertexindex_to_opposing_edgeindex( int v )
{
    if     ( v == 0 ) return 2;
    else if( v == 1 ) return 1;
    else if( v == 2 ) return 0;
    else assert(false);
}

int ManifoldTriangulation2D::edgeindex_to_opposing_vertexindex( int e )
{
    if     ( e == 0 ) return 2;
    else if( e == 1 ) return 1;
    else if( e == 2 ) return 0;
    else assert(false);
}

bool ManifoldTriangulation2D::vertexlists_equivalent( std::array<int,2> e1, std::array<int,2> e2 )
{
    assert( e1[0] != e1[1] && e2[0] != e2[1] );
    
    if     ( e1[0] == e2[0] && e1[1] == e2[1] ) return true;
    else if( e1[0] == e2[1] && e1[1] == e2[0] ) return true;
    else                                        return false;
}








// // int ManifoldTriangulation2D::is_vertexparents_cyclic( int v );
// // int ManifoldTriangulation2D::count_vertexparents( int v, int* ts );
// // int ManifoldTriangulation2D::list_vertexparents( int v, int* ts );
// // int ManifoldTriangulation2D::count_vertexparents_edges( int v, int* es );
// // int ManifoldTriangulation2D::list_vertexparents_edges( int v, int* es );



/* bisect a single edge */

void ManifoldTriangulation2D::bisect_edge( int e )
{
    assert( 0 <= e && e < counter_edges );
    
    if( count_edge_parents( e ) == 1 )
        bisect_outer_edge( e );
    else if( count_edge_parents( e ) == 2 )
        bisect_inner_edge( e );
    else
        assert(false);
}

void ManifoldTriangulation2D::bisect_inner_edge( int e )
{
    assert( 0 <= e && e < counter_edges );
    
    
}

void ManifoldTriangulation2D::bisect_outer_edge( int e )
{
    assert( 0 <= e && e < counter_edges );
    assert( count_edge_parents( e ) == 1 );
    
    int o = 1;
    
    FloatVector midpoint = get_edge_midpoint( e );
    coordinates.append( midpoint );
    
    int t       = get_edge_parents(e)[0];
    int el      = indexof_triangle_edge( t, e );
    int el_prev = get_prev_edge( el, o );
    int el_next = get_next_edge( el, o );
    
    assert( el != el_prev && el != el_next && el_next != el_prev );
    
    assert( data_triangle_edges[t][el] == get_triangle_edges(t)[el] );
    assert( e == get_triangle_edges(t)[el] );
    
    int e_prev = get_triangle_edges(t)[ el_prev ];
    int n_prev = get_triangle_neighbors(t)[ edgeindex_to_neighborindex(el_prev) ];
    int e_next = get_triangle_edges(t)[ el_next ];
    int n_next = get_triangle_neighbors(t)[ edgeindex_to_neighborindex(el_next) ];
    
    assert( e != e_prev && e != e_next && e_next != e_prev );
    assert( t != n_prev && t != n_next );
    if( n_prev != nullindex ) assert( n_next != n_prev );
    
    int vl_opp  = edgeindex_to_opposing_vertexindex( el );
    int vl_forw = edgeindex_to_opposing_vertexindex( el_prev );
    int vl_back = edgeindex_to_opposing_vertexindex( el_next );
    
    assert( vl_opp != vl_forw && vl_opp != vl_back && vl_forw != vl_back );
    
    std::array<int,3> t_new_vertices_back 
      = { get_triangle_vertices(t)[vl_opp ],
          get_triangle_vertices(t)[vl_back],
          counter_vertices
        };
    
    std::array<int,3> t_new_vertices_forward
      = { counter_vertices,
          get_triangle_vertices(t)[vl_forw],
          get_triangle_vertices(t)[vl_opp ]
        };
    
    std::array<int,3> t_new_edges_back;
    t_new_edges_back[ vertexindex_to_opposing_edgeindex(0) ] = e;
    t_new_edges_back[ vertexindex_to_opposing_edgeindex(1) ] = counter_edges + 1;
    t_new_edges_back[ vertexindex_to_opposing_edgeindex(2) ] = e_prev;
    
    std::array<int,3> t_new_edges_forward;
    t_new_edges_forward[ vertexindex_to_opposing_edgeindex(0) ] = e_next;
    t_new_edges_forward[ vertexindex_to_opposing_edgeindex(1) ] = counter_edges + 1;
    t_new_edges_forward[ vertexindex_to_opposing_edgeindex(2) ] = counter_edges;
    
    std::array<int,2> e_new_parents_back    = {                 t, nullindex };
    std::array<int,2> e_new_parents_forward = { counter_triangles, nullindex };
    std::array<int,2> e_new_pair            = { counter_triangles, t };
    
    
    data_triangle_vertices[t] = t_new_vertices_back;
    data_triangle_edges   [t] = t_new_edges_back;
    data_triangle_vertices.push_back( t_new_vertices_forward );
    data_triangle_edges.push_back   ( t_new_edges_forward );
    
    data_edge_parents[e] = e_new_parents_back;
    data_edge_parents.push_back( e_new_parents_forward );
    data_edge_parents.push_back( e_new_pair );
    
//     std::cout << "test: " << el << space << el_prev << space << el_next << nl;
//     std::cout << "test: " << e << space << e_prev << space << e_next << nl;
//     std::cout << "test: " << t << space << n_prev << space << n_next << nl;
    data_edge_parents[e_prev] = {                 t, n_prev };
    data_edge_parents[e_next] = { counter_triangles, n_next };
    
    counter_triangles += 1;
    counter_edges     += 2;
    counter_vertices  += 1;
    
}















FloatVector ManifoldTriangulation2D::get_triangle_midpoint( int t )
{
    assert( 0 <= t && t < counter_triangles );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = (   coordinates.getdata( get_triangle_vertices(t)[0], d )
                 + coordinates.getdata( get_triangle_vertices(t)[1], d )
                 + coordinates.getdata( get_triangle_vertices(t)[2], d )
               ) / 3.;
    return mid; 
}

FloatVector ManifoldTriangulation2D::get_edge_midpoint    ( int e )
{
    assert( 0 <= e && e < counter_edges );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = ( coordinates.getdata( get_edge_vertices(e)[0], d ) 
                 + coordinates.getdata( get_edge_vertices(e)[1], d )
               ) / 2.;
    return mid;
}


        



