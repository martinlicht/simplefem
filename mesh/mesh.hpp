#ifndef INCLUDEGUARD_MESH_ABSTRACTMESH
#define INCLUDEGUARD_MESH_ABSTRACTMESH


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
#include "../dense/functions.hpp"
#include "coordinates.hpp"


/*******************
****  
****  
****  Mesh Class 
****  
****  - contains Coordinate Class 
****  - intrinsic and exterior dimension
****  - provides relation about sub and supersimplices
****  - the topological data reflect a simplicial complex  
****    
****  
*******************/


class Mesh
{
    
    public:
        
        
        Mesh( int inner, int outer );
        
        void check() const;
        
        virtual void print( std::ostream& out ) const = 0;
        
        static const int nullindex; 
        
        static int is_not_nullindex( int i ){ return i != nullindex; }
        
        static int is_nullindex( int i ){ return i == nullindex; }
        
        
        
        /* Basic data */
        
        int getinnerdimension() const;
        
        int getouterdimension() const;
        
        Coordinates& getcoordinates();
        
        const Coordinates& getcoordinates() const;
        
        
        
        /* Static auxiliary functions */
        
        int index_from_pair( int sup, int sub ) const;
        
        void index_to_pair( int index, int& sup, int& sub ) const;
        
        int count_subsimplices( int sup, int sub ) const;
        
        
        
        /* Counting simplices */
        
        virtual bool dimension_counted( int dim ) const = 0;
        
        virtual int count_simplices( int dim ) const = 0;
        
        
        
        /* 
         * Accessing subsimplices 
         *
         * - check whether subsimplices are listed at all 
         * - get the subsimplex list of a cell
         * - test for subsimplex relation 
         * - get local index of subsimplex 
         * - get listed subsimplex 
         * 
         */
        
        virtual bool subsimplices_listed( int sup, int sub ) const = 0;
        
        virtual IndexMap getsubsimplices( int sup, int sub, int cell ) const = 0;
        
        virtual bool is_subsimplex( int sup, int sub, int cellsup, int cellsub ) const;
        
        virtual int get_subsimplex_index( int sup, int sub, int cellsup, int cellsub ) const;
        
        virtual int get_subsimplex( int sup, int sub, int cellsup, int localindex ) const;
        
        
        
        /* 
         * Accessing supersimplices 
         * 
         * - check whether supersimplices are listed at all 
         * - get the supersimplex list of a cell
         * - test for supersimplex relation 
         * - get local index of supersimplex 
         * 
         */
        
        virtual bool supersimplices_listed( int sup, int sub ) const = 0;
        
        virtual const std::vector<int> getsupersimplices( int sup, int sub, int cell ) const = 0;
        
        virtual bool is_supersimplex( int sup, int sub, int cellsup, int cellsub ) const;
        
        virtual int get_firstparent_of_subsimplex( int sup, int sub, int cellsub ) const;
        
        virtual int get_nextparent_of_subsimplex( int sup, int sub, int cellsup, int cellsub ) const;
        
        virtual int get_nextparent_by_localindex( int sup, int sub, int cellsup, int localindex ) const;
        
//         virtual IndexMap getnextparents( int sup, int sub, int cell ) const;
        
        virtual int get_index_of_supersimplex( int sup, int sub, int cellsup, int cellsub ) const;
        
        virtual int get_supersimplex_by_index( int sup, int sub, int cellsub, int parentindex ) const;
        
        // TODO: Iterator interface
        // ContainerInterface -- container.hpp 
        // derive from that class 
        
        
        /* 
         * Accessing geometric information about the mesh
         * 
         */
        
        Float getDiameter( int dim, int index ) const;
        
        Float getMeasure( int dim, int index ) const;
        
        Float getShapemeasure( int dim, int index ) const;
        Float getShapemeasure( int dim ) const;
        Float getShapemeasure() const;
        
        DenseMatrix getVertexCoordinateMatrix( int dim, int index ) const;
        
        DenseMatrix getTransformationJacobian( int dim, int index ) const;
        
        DenseMatrix getGradientProductMatrix( int dim, int index ) const;
        
        DenseMatrix getGradientProductMatrixRightFactor( int dim, int index ) const;
        
        
    private:
        
        int innerdimension;
        int outerdimension;
        
        Coordinates coordinates;
        
    private:
        
        std::map< std::pair<int,int>, std::vector<IndexMap> > auxdata;
        
    public: 
        
        const auto& getauxdata() { return auxdata; };
        
};




inline std::ostream& operator<<( std::ostream& os, const Mesh& mesh )
{
    mesh.print( os );
    return os;
}




// static inline int countsubsimplices( int n, int k )
// {
//     assert( 0 <= k && k <= n );
//     return binomial<int>( n+1, k+1 );
// }











#endif
