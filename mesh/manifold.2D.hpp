#ifndef INCLUDEGUARD_MANIFOLD_2D
#define INCLUDEGUARD_MANIFOLD_2D


#include <string>
#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>


#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"
#include "coordinates.hpp"


/*******************
****  
****  
****  ManifoldTriangulation2D Class 
****  
****  - specialized mesh class 
****  
****    Rationale:
****    - for every triangle, we save the 3 vertices 
****    - for every triangle, we save the 3 edges 
****    - for every edge, we save the 2 parents in ascending order 
****    - for every vertex, we save the first parent 
****     
****    Traversing the triangles of a vertex is tricky.
****    - if the enumeration in each vertex is oriented, 
****      then no additional data is needed, except the first triangle 
****    - alternatively:
****      we save the next triangle for each triangle 3x
****    - a good compromise:
****      save the first two triangles,
****      which dictates the transversal
****    From the definition of the Euler characteristic we deduce 
****    that in a typical triangulation we have got as many edges 
****    as we have got triangles and vertices. Moreover, 
****    bisecting an interior edge introduces more triangles 
****    than vertices, so we expect this require less memory.  
****  
*******************/


class ManifoldTriangulation2D
{

    public:
    
        ManifoldTriangulation2D( int outerdim );
        
        ManifoldTriangulation2D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,3>> triangles
        );
        
        virtual ~ManifoldTriangulation2D();
        
        virtual void check() const;

        virtual void print( std::ostream& out ) const;

        static const int nullindex = std::numeric_limits<int>::max();
        
        
        /* Coordinates */
        
        int getouterdimension() const;
        Coordinates& getcoordinates();
        const Coordinates& getcoordinates() const;
        
        /* General management */
        
        /* count the simplices of a certain type */
        int count_triangles() const;
        int count_edges() const;
        int count_vertices() const;
        
        /* contains a subsimplex? */
        bool contains_triangle_edge  ( int t, int e ) const; 
        bool contains_triangle_vertex( int t, int v ) const;
        bool contains_edge_vertex    ( int e, int v ) const;
        
        /* get index of a subsimplex */
        int indexof_triangle_edge  ( int t, int e ) const; 
        int indexof_triangle_vertex( int t, int v ) const;
        int indexof_edge_vertex    ( int e, int v ) const;
        
        /* get the subsimplices */
        const  std::array<int,3>  get_triangle_edges   ( int t ) const; 
        const  std::array<int,3>  get_triangle_vertices( int t ) const;
        const  std::array<int,2>  get_edge_vertices    ( int e ) const;
        
        
        /* Triangle neighbors */
        
        bool is_triangle_neighbor( int t, int nt ) const; 
        
        int indexof_triangle_neighbor( int t, int nt ) const; 
        
        const  std::array<int,3>  get_triangle_neighbors( int t ) const; 
        
        bool is_edge_between( int t1, int t2, int e ) const;
        
        int get_edge_between( int t1, int t2 ) const;
        
        
        
        /* parents of an edge */
        
        bool is_edge_parent( int t, int e ) const;
        
        int indexof_edge_parent( int t, int e ) const;
        
        int count_edge_parents( int e ) const;
        
        const  std::array<int,2>  get_edge_parents( int e ) const;
        
        
        
        
        /* index conversion */
        
        static int neighborindex_to_edgeindex( int );
        
        static int edgeindex_to_neighborindex( int );
        
        static std::array<int,2> edgeindex_to_vertexindices( int );
        
        static int vertexindex_to_opposing_edgeindex( int );
        
        static int edgeindex_to_opposing_vertexindex( int );
        
        
        
        
        
        
// //         int is_vertexparents_cyclic( int v );
// //         int count_vertexparents( int v, int* ts );
// //         int list_vertexparents( int v, int* ts );
// //         int count_vertexparents_edges( int v, int* es );
// //         int list_vertexparents_edges( int v, int* es );
        
        
        
        /* bisect a single edge */
        
        void bisect_edge( int e );
        
        void bisect_inner_edge( int e );
        
        void bisect_outer_edge( int e );
        
        
        
    private:

        int outerdimension;

        Coordinates coordinates;

        int counter_triangles;
        int counter_edges;
        int counter_vertices;
        
        std::vector< std::array<int,3> > data_triangle_vertices;
        std::vector< std::array<int,3> > data_triangle_edges;
        std::vector< std::array<int,2> > data_edge_parents;
//         std::vector< std::array<int,2> > data_vertex_firstparents;
        
        
};



static inline int countsubsimplices( int n, int k )
{
    assert( 0 <= k && k <= n );
    return binomial<int>( n+1, k+1 );
}




inline std::ostream& operator<<( std::ostream& os, const ManifoldTriangulation2D& mt2d )
{
    mt2d.print( os );
    return os;
}


#endif