#ifndef INCLUDEGUARD_MESH_SIMPLICIAL_ND_HPP
#define INCLUDEGUARD_MESH_SIMPLICIAL_ND_HPP


#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>


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
            const std::vector<int>& volume_vertices
        );
        
        explicit MeshSimplicialND( const Mesh& mesh );
        
        MeshSimplicialND( const MeshSimplicialND& mesh ) = default;
        
        virtual ~MeshSimplicialND();
        
        bool compare( const MeshSimplicialND& ) const;
        
        virtual void check() const;
        
        virtual void print( std::ostream& out ) const override;
        
        
        /* inherited methods */
        
        virtual bool dimension_counted( int dim ) const override;
        
        virtual int count_simplices( int dim ) const override;
        
        virtual bool subsimplices_listed( int sup, int sub ) const override;
        
        virtual IndexMap getsubsimplices( int sup, int sub, int cell ) const override;
        
        virtual bool supersimplices_listed( int sup, int sub ) const override;
        
        virtual const std::vector<int> getsupersimplices( int sup, int sub, int cell ) const override;
        
        
        virtual SimplexFlag get_flag( int dim, int index ) const override;
        
        virtual void set_flag( int dim, int index, SimplexFlag flag ) override;
        
        
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
        
        std::vector<std::vector<int>> data_subsimplices;
        std::vector<std::vector<int>> data_firstparents;
        std::vector<std::vector<int>> data_nextparents;
    
        std::vector<std::vector<SimplexFlag>> flags_simplices;
    
};




// inline std::ostream& operator<<( std::ostream& os, const MeshSimplicialND& mND )
// {
//     mND.print( os );
//     return os;
// }



inline bool operator==( const MeshSimplicialND& m1, const MeshSimplicialND& m2 )
{
    return m1.compare( m2 );
}

inline bool operator!=( const MeshSimplicialND& m1, const MeshSimplicialND& m2 )
{
    return !( m1 == m2 );
}






#endif
