#ifndef INCLUDEGUARD_MESH_SIMPLICIAL_ND
#define INCLUDEGUARD_MESH_SIMPLICIAL_ND


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
****  MeshSimplicialND Class 
****  
****    
****    
*******************/


class MeshSimplicialND
: public Mesh
{

    public:
    
        MeshSimplicialND( int innerdim, int outerdim );
        
        MeshSimplicialND( 
            int innerdim,
            int outerdim,
            const Coordinates& coords,
            const std::vector<int> volume_vertices
        );
        
        virtual ~MeshSimplicialND();
        
        bool operator== ( const MeshSimplicialND& ) const;
        
        bool operator!= ( const MeshSimplicialND& ) const;
        
        virtual void check() const;
        
        virtual void print( std::ostream& out ) const override;
        
        
        /* inherited methods */
        
        virtual bool dimensioncounted( int dim ) const override;
        
        virtual int countsimplices( int dim ) const override;
        
        virtual bool subsimplices_listed( int sup, int sub ) const override;
        
        virtual const IndexMap getsubsimplices( int sup, int sub, int cell ) const override;
        
        virtual bool supersimplices_listed( int sup, int sub ) const override;
        
        virtual const std::vector<int> getsupersimplices( int sup, int sub, int cell ) const override;
        
        
        /* refinement */
        
        void bisect_edge( int e );
        
        void uniformrefinement();
        
        /* other things */
        
        FloatVector get_simplex_midpoint( int dim, int cell ) const;
        
        
    private: 
        
        void rebuild();
        
        int relationindex( int sup, int sub ) const;
        
        int subsimplexnumber( int sup, int sub ) const;
        
        
    private:
        
        std::vector<int> counter_simplices;
        
        std::vector<std::vector<int>> data_simplex_subsimplices;
        std::vector<std::vector<int>> data_simplex_firstparent_simplex;
        std::vector<std::vector<int>> data_simplex_nextparent_of_subsimplices;
    
    
};




inline std::ostream& operator<<( std::ostream& os, const MeshSimplicialND& mND )
{
    mND.print( os );
    return os;
}




inline MeshSimplicialND UnitSquare()
{
    return MeshSimplicialND(
      1, 3, 
      Coordinates( 3, 2, {
        -1., 0., 1., // 0
         1., 3., 2., // 1
      } ),
      {
        0, 1 
      }
    );
}

inline MeshSimplicialND StandardSquare()
{
    return MeshSimplicialND(
      1, 1, 
      Coordinates( 1, 2, {
         0., // 0
         1., // 1
      } ),
      {
        0, 1
      }
    );
}











#endif
