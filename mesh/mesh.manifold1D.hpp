#ifndef INCLUDEGUARD_MESH_MANIFOLD_1D
#define INCLUDEGUARD_MESH_MANIFOLD_1D


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
****  MeshManifold1D Class 
****  
****  - specialized mesh class for one-dimensional meshes 
****  - looks like an undirected graph. Parallel edges are allowed.
****  - We want to get from every edge to its vertices 
****    and from every vertex we want to traverse its parent edges 
****  
****    Rationale:
****    - for every edge, we save the 2 vertices 
****    - for every vertex, we save the first parent edge
****    - for every edge, we save the next parents of each vertex. 
****    
****    
*******************/


class MeshManifold1D
: public Mesh
{

    public:
    
        MeshManifold1D( int outerdim = 1 );
        
        MeshManifold1D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,2>> edge_vertices
        );
        
        MeshManifold1D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,2>> edge_vertices,
            const std::vector<int              > vertex_firstparent,
            const std::vector<std::array<int,2>> edge_nextparents
        );
        
        
        virtual ~MeshManifold1D();
        
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
        int count_edges()    const;
        int count_vertices() const;
        
        /* subsimplex relation of edges and vertices */
        
        bool contains_edge_vertex( int e, int v ) const;
        
        int indexof_edge_vertex( int e, int v ) const;
        
        const std::array<int,2> get_edge_vertices ( int e ) const;
        
        
        /* edge parents of a vertex */
        
        int count_vertex_edge_parents( int v ) const;
        
        int get_vertex_firstparent( int v ) const;
        
        int get_vertex_nextparent( int v, int e ) const;
        
        bool is_edge_vertex_parent( int e, int v ) const;
        
        int indexof_edge_vertex_parent( int e, int v ) const;
        
        std::vector<int> get_edge_parents_of_vertex( int v ) const;
        
        
        /* refinement */
        
        void refineedge( int e );
        
        void uniformrefinement();

        
        /* other things */
        
        FloatVector get_edge_midpoint( int e );
        
        
    private:

        int counter_edges;
        int counter_vertices;
        
        std::vector< std::array<int,2> > data_edge_vertices;
        std::vector< int               > data_vertex_firstparent;
        std::vector< std::array<int,2> > data_edge_nextparents;
        
        
        
};




inline std::ostream& operator<<( std::ostream& os, const MeshManifold1D& mt2d )
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




inline MeshManifold1D UnitSquare()
{
    return MeshManifold1D(
      1,
      Coordinates( 1, 2, {
        -1., // 0
         1., // 1
      } ),
      {
        { 0, 1 }
      }
    );
}

inline MeshManifold1D StandardSquare()
{
    return MeshManifold1D(
      1,
      Coordinates( 1, 2, {
        -1., // 0
         1., // 1
      } ),
      {
        { 0, 1 }
      }
    );
}











#endif
