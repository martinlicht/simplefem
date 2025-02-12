#ifndef INCLUDEGUARD_MESH_SIMPLICIAL_ND_HPP
#define INCLUDEGUARD_MESH_SIMPLICIAL_ND_HPP


// #include <ostream>
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
    
        /* Constructors */
        
        MeshSimplicialND( int innerdim, int outerdim );
        
        MeshSimplicialND( 
            int innerdim,
            int outerdim,
            const Coordinates& coords,
            const std::vector<int>& volume_vertices
        );
        
        explicit MeshSimplicialND( const Mesh& mesh );
        
        /* standard interface */
        
        MeshSimplicialND( const MeshSimplicialND& mesh ) = default;
        
        virtual ~MeshSimplicialND();
        
        /* standard methods for operators */
        
        virtual void check() const;
        
        // virtual void print( std::ostream& out ) const override;

        virtual std::string text() const override ;
        
        /* OTHER METHODS */
        
        bool is_equal_to( const MeshSimplicialND& ) const;
        
        
        /* inherited methods */
        
        virtual bool has_dimension_counted( int dim ) const override;
        
        virtual int count_simplices( int dim ) const override;
        
        virtual bool has_subsimplices_listed( int sup, int sub ) const override;
        
        virtual IndexMap getsubsimplices( int sup, int sub, int cell ) const override;
        
        virtual bool has_supersimplices_listed( int sup, int sub ) const override;
        
        virtual const std::vector<int> getsupersimplices( int sup, int sub, int cell ) const override;
        
        
        virtual SimplexFlag get_flag( int dim, int index ) const override;
        
        virtual void set_flag( int dim, int index, SimplexFlag flag ) override;
        
        
        /* refinement */
        
        void bisect_edge( int e );
        
        void uniformrefinement();
        
        /* other things */
        
        FloatVector get_simplex_midpoint( int dim, int cell ) const;

        std::size_t memorysize() const override;
        
        
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


    public:

        inline bool operator==( const MeshSimplicialND& m2 ) const 
        {
            return this->is_equal_to( m2 );
        }

        inline bool operator!=( const MeshSimplicialND& m2 ) const 
        {
            return !( *this == m2 );
        }

    
};













#endif
