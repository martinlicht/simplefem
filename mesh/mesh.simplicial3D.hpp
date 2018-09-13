#ifndef INCLUDEGUARD_MESH_SIMPLICIAL_3D
#define INCLUDEGUARD_MESH_SIMPLICIAL_3D


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
****  MeshSimplicial3D Class 
****  
****  - specialized mesh class for finite two-dimensional simplicial complexes
****    
****    
*******************/


class MeshSimplicial3D
: public Mesh
{

    public:
    
        explicit MeshSimplicial3D( int outerdim = 3 );
        
        MeshSimplicial3D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,4>>& tetrahedron_vertices
        );
        
        MeshSimplicial3D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,4>>& tetrahedron_faces,
            const std::vector<int              >& face_firstparent_tetrahedron,
            const std::vector<std::array<int,4>>& tetrahedron_nextparents_of_faces,
            const std::vector<std::array<int,6>>& tetrahedron_edges,
            const std::vector<int              >& edge_firstparent_tetrahedron,
            const std::vector<std::array<int,6>>& tetrahedron_nextparents_of_edges,
            const std::vector<std::array<int,4>>& tetrahedron_vertices,
            const std::vector<int              >& vertex_firstparent_tetrahedron,
            const std::vector<std::array<int,4>>& tetrahedron_nextparents_of_vertices,
            const std::vector<std::array<int,3>>& face_edges,
            const std::vector<int              >& edge_firstparent_face,
            const std::vector<std::array<int,3>>& face_nextparents_of_edges,
            const std::vector<std::array<int,3>>& face_vertices,
            const std::vector<int              >& vertex_firstparent_face,
            const std::vector<std::array<int,3>>& face_nextparents_of_vertices,
            const std::vector<std::array<int,2>>& edge_vertices,
            const std::vector<int              >& vertex_firstparent_edge,
            const std::vector<std::array<int,2>>& edge_nextparents_of_vertices
        );
        
        
        virtual ~MeshSimplicial3D();
        
        virtual void check() const;
        
        virtual void print( std::ostream& out ) const override;
        
        bool operator== ( const MeshSimplicial3D& mesh ) const; 
        
        bool operator!= ( const MeshSimplicial3D& mesh ) const; 
        
        /* inherited methods */
        
        virtual bool dimension_counted( int dim ) const override;
        
        virtual int count_simplices( int dim ) const override;
        
        virtual bool subsimplices_listed( int sup, int sub ) const override;
        
        virtual IndexMap getsubsimplices( int sup, int sub, int cell ) const override;
        
        virtual bool supersimplices_listed( int sup, int sub ) const override;
        
        virtual const std::vector<int> getsupersimplices( int sup, int sub, int cell ) const override;
        
        
        /* General management */
        
        
        /* count the simplices of a certain type */
        int count_tetrahedra() const;
        int count_faces()      const;
        int count_edges()      const;
        int count_vertices()   const;
        
        
        /* subsimplex relation of tetrahedra and faces */
        
        bool contains_tetrahedron_face( int t, int f ) const;
        
        int indexof_tetrahedron_face( int t, int f ) const;
        
        int get_tetrahedron_face( int t, int fi ) const;
        
        const std::array<int,4> get_tetrahedron_faces( int t ) const;
        
        
        /* subsimplex relation of tetrahedra and edges */
        
        bool contains_tetrahedron_edge( int t, int e ) const;
        
        int indexof_tetrahedron_edge( int t, int e ) const;
        
        int get_tetrahedron_edge( int t, int ei ) const;
        
        const std::array<int,6> get_tetrahedron_edges( int t ) const;
        
        
        /* subsimplex relation of tetrahedra and vertices */
        
        bool contains_tetrahedron_vertex( int t, int v ) const;
        
        int indexof_tetrahedron_vertex( int t, int v ) const;
        
        int get_tetrahedron_vertex( int t, int vi ) const;
        
        const std::array<int,4> get_tetrahedron_vertices ( int t ) const;
        
        
        
        /* subsimplex relation of faces and edges */
        
        bool contains_face_edge( int f, int e ) const;
        
        int indexof_face_edge( int f, int e ) const;
        
        int get_face_edge( int f, int ei ) const;
        
        const std::array<int,3> get_face_edges( int f ) const;
        
        
        /* subsimplex relation of faces and vertices */
        
        bool contains_face_vertex( int f, int v ) const;
        
        int indexof_face_vertex( int f, int v ) const;
        
        int get_face_vertex( int f, int vi ) const;
        
        const std::array<int,3> get_face_vertices ( int f ) const;
        
        
        /* subsimplex relation of edges and vertices */
        
        bool contains_edge_vertex( int e, int v ) const;
        
        int indexof_edge_vertex( int e, int v ) const;
        
        int get_edge_vertex( int e, int vi ) const;
        
        const std::array<int,2> get_edge_vertices ( int e ) const;
        
        
        
        
        
        
        
        
        
        /* tetrahedron parents of a face */
        
        int count_face_tetrahedron_parents( int f ) const;
        
        int get_face_firstparent_tetrahedron( int f ) const;
        
        int get_face_nextparent_tetrahedron( int f, int t ) const;
        
        int get_tetrahedron_nextparent_of_face( int t, int fi ) const;
        
        bool is_tetrahedron_face_parent( int t, int f ) const;
        
        int indexof_tetrahedron_face_parent( int t, int f ) const;
        
        std::vector<int> get_tetrahedron_parents_of_face( int f ) const;
        
        
        /* tetrahedron parents of an edge */
        
        int count_edge_tetrahedron_parents( int e ) const;
        
        int get_edge_firstparent_tetrahedron( int e ) const;
        
        int get_edge_nextparent_tetrahedron( int e, int t ) const;
        
        int get_tetrahedron_nextparent_of_edge( int t, int ei ) const;
        
        bool is_tetrahedron_edge_parent( int t, int e ) const;
        
        int indexof_tetrahedron_edge_parent( int t, int e ) const;
        
        std::vector<int> get_tetrahedron_parents_of_edge( int e ) const;
        
        
        /* tetrahedron parents of a vertex */
        
        int count_vertex_tetrahedron_parents( int v ) const;
        
        int get_vertex_firstparent_tetrahedron( int v ) const;
        
        int get_vertex_nextparent_tetrahedron( int v, int t ) const;
        
        int get_tetrahedron_nextparent_of_vertex( int t, int vi ) const;
        
        bool is_tetrahedron_vertex_parent( int t, int v ) const;
        
        int indexof_tetrahedron_vertex_parent( int t, int v ) const;
        
        std::vector<int> get_tetrahedron_parents_of_vertex( int v ) const;
        
        
        /* face parents of an edge */
        
        int count_edge_face_parents( int e ) const;
        
        int get_edge_firstparent_face( int e ) const;
        
        int get_edge_nextparent_face( int e, int f ) const;
        
        int get_face_nextparent_of_edge( int f, int ei ) const;
        
        bool is_face_edge_parent( int f, int e ) const;
        
        int indexof_face_edge_parent( int f, int e ) const;
        
        std::vector<int> get_face_parents_of_edge( int e ) const;
        
        
        /* face parents of a vertex */
        
        int count_vertex_face_parents( int v ) const;
        
        int get_vertex_firstparent_face( int v ) const;
        
        int get_vertex_nextparent_face( int v, int f ) const;
        
        int get_face_nextparent_of_vertex( int f, int vi ) const;
        
        bool is_face_vertex_parent( int f, int v ) const;
        
        int indexof_face_vertex_parent( int f, int v ) const;
        
        std::vector<int> get_face_parents_of_vertex( int v ) const;
        
        
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
        
        void longest_edge_bisection( std::vector<int> edges );
        
        void uniformrefinement();
        
        void midpoint_refinement( int t );
        
        void midpoint_refinement_global();
        
        
        /* other things */
        
        FloatVector get_tetrahedron_midpoint( int t ) const;
        FloatVector get_face_midpoint( int f ) const;
        FloatVector get_edge_midpoint( int e ) const;
        Float get_edge_length( int e ) const;
        
        
    private:

        int counter_tetrahedra;
        int counter_faces;
        int counter_edges;
        int counter_vertices;
        
        std::vector< std::array<int,4> > data_tetrahedron_faces;
        std::vector< int               > data_face_firstparent_tetrahedron;
        std::vector< std::array<int,4> > data_tetrahedron_nextparents_of_faces;
        
        std::vector< std::array<int,6> > data_tetrahedron_edges;
        std::vector< int               > data_edge_firstparent_tetrahedron;
        std::vector< std::array<int,6> > data_tetrahedron_nextparents_of_edges;
        
        std::vector< std::array<int,4> > data_tetrahedron_vertices;
        std::vector< int               > data_vertex_firstparent_tetrahedron;
        std::vector< std::array<int,4> > data_tetrahedron_nextparents_of_vertices;
        
        std::vector< std::array<int,3> > data_face_edges;
        std::vector< int               > data_edge_firstparent_face;
        std::vector< std::array<int,3> > data_face_nextparents_of_edges;
        
        std::vector< std::array<int,3> > data_face_vertices;
        std::vector< int               > data_vertex_firstparent_face;
        std::vector< std::array<int,3> > data_face_nextparents_of_vertices;
        
        std::vector< std::array<int,2> > data_edge_vertices;
        std::vector< int               > data_vertex_firstparent_edge;
        std::vector< std::array<int,2> > data_edge_nextparents_of_vertices;

};




// inline std::ostream& operator<<( std::ostream& os, const MeshSimplicial3D& mt2d )
// {
//     mt2d.print( os );
//     return os;
// }








inline MeshSimplicial3D UnitSimplex3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         0., 0., 0., // 0
         0., 0., 1., // 1
         0., 1., 0., // 2
         1., 0., 0.  // 3
      } ),
      {
        { 0, 1, 2, 3 }
      }
    );
}


inline MeshSimplicial3D OctoDiamond3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         0., 0., 0., // 0
         1., 0., 0., // 1
         0., 1., 0., // 2
        -1., 0., 0., // 3
         0.,-1., 0., // 4
         0., 0., 1., // 5
         0., 0.,-1.  // 6
      } ),
      {
        { 0, 1, 2, 5 },
        { 0, 2, 3, 5 },
        { 0, 3, 4, 5 },
        { 0, 1, 4, 5 },
        { 0, 1, 2, 6 },
        { 0, 2, 3, 6 },
        { 0, 3, 4, 6 },
        { 0, 1, 4, 6 } 
      }
    );
}


inline MeshSimplicial3D StandardCube3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 8, {
         0., 0., 0., // 0
         0., 0., 1., // 0
         0., 1., 0., // 0
         0., 1., 1., // 0
         1., 0., 0., // 0
         1., 0., 1., // 0
         1., 1., 0., // 0
         1., 1., 1.  // 0
      } ),
      {
        { 0b000, 0b100, 0b110, 0b111 }, 
        { 0b000, 0b100, 0b101, 0b111 },
        { 0b000, 0b010, 0b110, 0b111 },
        { 0b000, 0b010, 0b011, 0b111 },
        { 0b000, 0b001, 0b101, 0b111 },
        { 0b000, 0b001, 0b011, 0b111 }
      }
    );
}








#endif
