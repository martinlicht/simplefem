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
        
        #ifdef NDEBUG
        static const int nullindex = std::numeric_limits<int>::max();
        #else
        static const int nullindex = 777777; 
        #endif
        
        static int is_not_nullindex( int i ){ return i != nullindex; }
        
        static int is_nullindex( int i ){ return i == nullindex; }
        
        
        
        /* Basic data */
        
        int getinnerdimension() const;
        
        int getouterdimension() const;
        
        Coordinates& getcoordinates();
        
        const Coordinates& getcoordinates() const;
        
        
        
        /* Counting simplices */
        
        virtual bool dimensioncounted( int dim ) const = 0;
        
        virtual int countsimplices( int dim ) const = 0;
        
        
        
        /* 
         * Accessing subsimplices 
         *
         * - check whether subsimplices are listed at all 
         * - get the subsimplex list of a cell
         * - test for subsimplex relation 
         * - get local index of subsimplex 
         * 
         */
        
        virtual bool subsimplices_listed( int sup, int sub ) const = 0;
        
        virtual const IndexMap getsubsimplices( int sup, int sub, int cell ) const = 0;
        
        bool is_subsimplex( int sup, int sub, int cellsup, int cellsub ) const;
        
        int get_subsimplix_index( int sup, int sub, int cellsup, int cellsub ) const;
        
        
        
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
        
        virtual int get_supersimplex_index( int sup, int sub, int cellsup, int cellsub ) const;
        
        // TODO: Iterator interface
        // ContainerInterface -- container.hpp 
        // derive from that class 
        
        
        
        
    private:
        
        int innerdimension;
        int outerdimension;
        
        Coordinates coordinates;
        
    private:
        
        std::map< std::pair<int,int>, std::vector<IndexMap> > auxdata;
        
    public: 
        
        const auto& getauxdata() { return auxdata; };
        
};






static inline int countsubsimplices( int n, int k )
{
    assert( 0 <= k && k <= n );
    return binomial<int>( n+1, k+1 );
}












inline std::ostream& operator<<( std::ostream& os, const Mesh& sm )
{
    sm.print( os );
    return os;
}


#endif