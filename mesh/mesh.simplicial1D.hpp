#ifndef INCLUDEGUARD_MESH_SIMPLICIAL_1D
#define INCLUDEGUARD_MESH_SIMPLICIAL_1D


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
****  MeshSimplicial1D Class 
****  
****  - specialized mesh class for finite one-dimensional simplicial complexes 
****  - looks like an undirected graph without loops or parallel edges.
****  
****    Content:
****    - for every edge, we save the 2 vertices 
****    - for every vertex, we save the first parent edge
****    - for every edge, we save the next parents of each vertex. 
****    
****    
*******************/


class MeshSimplicial1D
: public Mesh
{

    public:
    
        MeshSimplicial1D( int outerdim = 1 );
        
        MeshSimplicial1D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,2>> edge_vertices
        );
        
        MeshSimplicial1D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,2>> edge_vertices,
            const std::vector<std::array<int,2>> edge_nextparents_of_vertices,
            const std::vector<int              > vertex_firstparent_edge
        );
        
        
        virtual ~MeshSimplicial1D();
        
        bool operator== ( const MeshSimplicial1D& ) const;
        
        bool operator!= ( const MeshSimplicial1D& ) const;
        
        virtual void check() const;
        
        virtual void print( std::ostream& out ) const override;
        
        
        /* inherited methods */
        
        virtual bool dimensioncounted( int dim ) const override;
        
        virtual int countsimplices( int dim ) const override;
        
        virtual bool subsimplices_listed( int sup, int sub ) const override;
        
        virtual IndexMap getsubsimplices( int sup, int sub, int cell ) const override;
        
        virtual bool supersimplices_listed( int sup, int sub ) const override;
        
        virtual const std::vector<int> getsupersimplices( int sup, int sub, int cell ) const override;
        
        
        /* General management */
        
        
        /* count the simplices of a certain type */
        int count_edges()    const;
        int count_vertices() const;
        
        
        /* subsimplex relation of edges and vertices */
        
        bool contains_edge_vertex( int e, int v ) const;
        
        int indexof_edge_vertex( int e, int v ) const;
        
        int get_edge_vertex( int e, int vi ) const;
        
        const std::array<int,2> get_edge_vertices ( int e ) const;
        
        
        /* edge parents of a vertex */
        
        int count_vertex_edge_parents( int v ) const;
        
        int get_vertex_firstparent_edge( int v ) const;
        
        int get_vertex_nextparent_edge( int v, int e ) const;
        
        int get_edge_nextparent_of_vertex( int e, int vi ) const;
        
        bool is_edge_vertex_parent( int e, int v ) const;
        
        int indexof_edge_vertex_parent( int e, int v ) const;
        
        std::vector<int> get_edge_parents_of_vertex( int v ) const;
        
        
        /* refinement */
        
        void bisect_edge( int e );
        
        void uniformrefinement();
        
        void improved_uniformrefinement();
        
        
        /* other things */
        
        FloatVector get_edge_midpoint( int e ) const;
        
        
    private:

        int counter_edges;
        int counter_vertices;
        
        std::vector< std::array<int,2> > data_edge_vertices;
        std::vector< int               > data_vertex_firstparent_edge;
        std::vector< std::array<int,2> > data_edge_nextparents_of_vertices;
        
        
        
};




// inline std::ostream& operator<<( std::ostream& os, const MeshSimplicial1D& mt1d )
// {
//     mt1d.print( os );
//     return os;
// }




inline MeshSimplicial1D UnitSquare()
{
    return MeshSimplicial1D(
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

inline MeshSimplicial1D StandardSquare()
{
    return MeshSimplicial1D(
      1,
      Coordinates( 1, 2, {
         0., // 0
         1., // 1
      } ),
      {
        { 0, 1 }
      }
    );
}











#endif
