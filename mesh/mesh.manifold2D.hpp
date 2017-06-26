#ifndef INCLUDEGUARD_MESH_MANIFOLD_2D
#define INCLUDEGUARD_MESH_MANIFOLD_2D


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
#include "mesh.hpp"


/*******************
****  
****  
****  MeshManifold2D Class 
****  
****  - specialized mesh class 
****  - 
****  
****  
****   
****    Rationale:
****    - for every triangle, we save the 3 vertices 
****    - for every triangle, we save the 3 edges 
****    - for every edge, we save the 2 parents in ascending order 
****    --- for every vertex, we save the first parent 
****     
****    Traversing the triangles of a vertex is tricky.
****    
****    - if the enumeration in each vertex is oriented, 
****      then no additional data is needed, except the first triangle 
****    - alternatively:
****      we save the next triangle for each triangle 3x
****    - a good compromise:
****      save the first two triangles, which dictates the transversal
****    - Finally, we can fix an orientation for each triangle,
****      which could e.g. be computed on-the-fly from the vertex enumeration,
****      and base our traversal on that. 
****    
****    From the definition of the Euler characteristic we deduce 
****    that in a typical triangulation we have got as many edges 
****    as we have got triangles and vertices. Moreover, 
****    bisecting an interior edge introduces more triangles 
****    than vertices, so we expect this require less memory.  
****    
****    
****    
****    
****    
*******************/


class MeshManifold2D
: public Mesh
{

    public:
    
        MeshManifold2D( int outerdim = 2 );
        
        MeshManifold2D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,3>> triangle_vertices
        );
        
        MeshManifold2D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,3>> triangle_vertices,
            const std::vector<std::array<int,3>> triangle_edges,
            const std::vector<std::array<int,2>> edge_parents
        );
        
        
        virtual ~MeshManifold2D();
        
        virtual void check() const;
        
        virtual void print( std::ostream& out ) const override;
        
        
        /* inherited methods */
        
        virtual bool dimensioncounted( int dim ) const override;
        
        virtual int countsimplices( int dim ) const override;
        
        virtual bool subsimplices_listed( int sup, int sub ) const override;
        
        virtual const IndexMap getsubsimplices( int sup, int sub, int cell ) const override;
        
        virtual bool supersimplices_listed( int sup, int sub ) const override;
        
        virtual const std::vector<int> getsupersimplices( int sup, int sub, int cell ) const override;
        
        
        
        
        /* General management */
        
        /* count the simplices of a certain type */
        int count_triangles() const;
        int count_edges()     const;
        int count_vertices()  const;
        
        /* contains a subsimplex? */
        bool contains_triangle_edge  ( int t, int e ) const; 
        bool contains_triangle_vertex( int t, int v ) const;
        bool contains_edge_vertex    ( int e, int v ) const;
        
        /* get index of a subsimplex */
        int indexof_triangle_edge  ( int t, int e ) const; 
        int indexof_triangle_vertex( int t, int v ) const;
        int indexof_edge_vertex    ( int e, int v ) const;
        
        /* get the subsimplices */
        const std::array<int,3> get_triangle_edges   ( int t ) const; 
        const std::array<int,3> get_triangle_vertices( int t ) const;
        const std::array<int,2> get_edge_vertices    ( int e ) const;
        
        
        /* Triangle neighbors */
        
        bool is_triangle_neighbor( int t, int nt ) const; 
        int indexof_triangle_neighbor( int t, int nt ) const; 
        const std::array<int,3> get_triangle_neighbors( int t ) const; 
        
        bool is_edge_between( int t1, int t2, int e ) const;
        int get_edge_between( int t1, int t2        ) const;
        
        
        /* orientation-related things */
        
        int get_forwarding_edge  ( int t, int el, int o );
        int get_backwarding_edge ( int t, int el, int o );
                
        static int get_prev_edge    ( int el, int o );
        static int get_next_edge    ( int el, int o );
        static int get_prev_neighbor( int nl, int o );
        static int get_next_neighbor( int nl, int o );
        static int get_prev_vertex  ( int vl, int o );
        static int get_next_vertex  ( int vl, int o );
        
        static int get_first_vertex( int o );
        static int get_last_vertex ( int o );
        
        
        
        /* triangle parents of an edge */
        
        bool is_triangle_edge_parent( int t, int e ) const;
        
        int indexof_triangle_edge_parent( int t, int e ) const;
        
        int count_triangle_edge_parents( int e ) const;
        
        const std::array<int,2> get_triangle_edge_parents( int e ) const;
        
        int get_edge_othertriangleparent( int t, int e ) const;
        
        int orientation_induced( int t, int el ) const; 
        
        
        
        
        /* triangle parents of a vertex */
        
        int indexof_triangle_vertex_parent( int t, int v ) const;
        
        int count_vertex_triangle_parents( int v ) const;
        
        const std::vector<int> get_triangle_parents_of_vertex( int v ) const;
        
        
        /* edge parents of a vertex */
        
        int indexof_edge_vertex_parent( int e, int v ) const;
        
        int count_vertex_edge_parents( int v ) const;
        
        const std::vector<int> get_edge_parents_of_vertex( int v ) const;
        
        
        
        
        
        /* index conversion */
        
        static int neighborindex_to_edgeindex( int n );
        
        static int edgeindex_to_neighborindex( int e );
        
        static std::array<int,2> edgeindex_to_vertexindices( int e );
        
        static std::array<int,2> duple_from_triple( std::array<int,3> triple, std::array<int,2> duple );
        
        static int vertexindex_to_opposing_edgeindex( int v );
        
        static int edgeindex_to_opposing_vertexindex( int e );
        
        static bool vertexpairs_equivalent( std::array<int,2> vp1, std::array<int,2> vp2 );
        
        
        
        
        /* uniform refinement */
        
        void uniformrefinement();
        
        
        
        /* bisect a single edge */
        
        void bisect_edge( int e );
        
        void bisect_inner_edge( int e );
        
        void bisect_outer_edge( int e );
        
        
        /* other things */
        
        FloatVector get_triangle_midpoint( int t );
        FloatVector get_edge_midpoint    ( int e );
        
        
        
        
    private:

        int counter_triangles;
        int counter_edges;
        int counter_vertices;
        
        std::vector< std::array<int,3> > data_triangle_vertices;
        std::vector< std::array<int,3> > data_triangle_edges;
        std::vector< std::array<int,2> > data_edge_parents;
        
        struct orientedtriangle{ int parenttriangle; signed char orientation; };
        
        std::vector< orientedtriangle > data_vertex_firstparents;
        
        
};




inline std::ostream& operator<<( std::ostream& os, const MeshManifold2D& mt2d )
{
    mt2d.print( os );
    return os;
}



inline bool duple_equivalent( std::array<int,2> d1, std::array<int,2> d2 )
{
  if( d1[0] == d2[0] && d1[1] == d2[1] ) return true;
  if( d1[0] == d2[1] && d1[1] == d2[0] ) return true;
  return false;
}

inline bool triple_equivalent( std::array<int,3> d1, std::array<int,3> d2 )
{
  if( d1[0] == d2[0] && d1[1] == d2[1] && d1[2] == d2[2] ) return true;
  if( d1[0] == d2[0] && d1[2] == d2[1] && d1[1] == d2[2] ) return true;
  if( d1[1] == d2[0] && d1[0] == d2[1] && d1[2] == d2[2] ) return true;
  if( d1[1] == d2[0] && d1[2] == d2[1] && d1[0] == d2[2] ) return true;
  if( d1[2] == d2[0] && d1[0] == d2[1] && d1[1] == d2[2] ) return true;
  if( d1[2] == d2[0] && d1[1] == d2[1] && d1[0] == d2[2] ) return true;
  return false;
}






inline MeshManifold2D UnitSquare()
{
    return MeshManifold2D(
      2,
      Coordinates( 2, 4, {
        -1., -1., // 0
        -1.,  1., // 1
         1., -1., // 2
         1.,  1.  // 3
      } ),
      {
        { 0, 1, 3 },
        { 0, 2, 3 }
      }
    );
}

inline MeshManifold2D StandardSquare()
{
    return MeshManifold2D(
      2,
      Coordinates( 2, 9, {
         0.,  0., // 0
         1.,  0., // 1
         1.,  1., // 2
         0.,  1., // 3
        -1.,  1., // 4
        -1.,  0., // 5
        -1., -1., // 6
         0., -1., // 7
         1., -1.  // 8
      } ),
      {
        { 0, 1, 2 },
        { 0, 2, 3 },
        { 0, 3, 4 },
        { 0, 4, 5 },
        { 0, 5, 6 },
        { 0, 6, 7 },
        { 0, 7, 8 },
        { 0, 8, 1 }
      }
    );
}



inline MeshManifold2D LShapedDomain()
{
    return MeshManifold2D(
      2,
      Coordinates( 2, 8, {
         0.,  0., // 0
         1.,  0., // 1
         1.,  1., // 2
         0.,  1., // 3
        -1.,  1., // 4
        -1.,  0., // 5
        -1., -1., // 6
         0., -1.  // 7
      } ),
      {
        { 0, 1, 2 },
        { 0, 2, 3 },
        { 0, 3, 4 },
        { 0, 4, 5 },
        { 0, 5, 6 },
        { 0, 6, 7 }
      }
    );
}



inline MeshManifold2D SlitDomain()
{
    return MeshManifold2D(
      3,
      Coordinates( 3, 10, {
         0.,  0., 0, // 0
         1.,  0., 0.5, // 1
         1.,  1., 0, // 2
         0.,  1., 0, // 3
        -1.,  1., 0, // 4
        -1.,  0., 0, // 5
        -1., -1., 0, // 6
         0., -1., 0, // 7
         1., -1., 0, // 8
         1.,  0., -0.5,  // 9
      } ),
      {
        { 0, 1, 2 }, // 0
        { 0, 2, 3 }, // 1
        { 0, 3, 4 }, // 2
        { 0, 4, 5 }, // 3
        { 0, 5, 6 }, // 4
        { 0, 6, 7 }, // 5
        { 0, 7, 8 }, // 6
        { 0, 8, 9 }  // 7
      }
    );
}









//         void get_firstvertexparent_triangle( int v, int& t, bool& orientation );
//         void get_nextvertexparent_triangle ( int v, int& t, bool& orientation );
//         void is_lastvertexparent_triangle  ( int v, int& t, bool& orientation );
//         
//         void get_firstvertexparent_edges( int v, int& e, bool& orientation ); // orientation is not enough
//         void get_nextvertexparent_edges ( int v, int& e, bool& orientation );
//         void is_lastvertexparent_edges  ( int v, int& e, bool& orientation );
//         
//         int are_vertexparents_cyclic     ( int v );
//         int count_vertexparents_triangles( int v );
//         int count_vertexparents_edges    ( int v );
//         
//         void initialize_vertexparents();
//         void initialize_vertexparents( int t );



#endif
