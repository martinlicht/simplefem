#ifndef INCLUDEGUARD_SIMPLICIALMESH
#define INCLUDEGUARD_SIMPLICIALMESH


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
****  SimplicialMesh Class 
****  
****  - contains Coordinate Class 
****  - intrinsic and exterior dimension
****  - provides relation about sub and supersimplices
****  
*******************/


class SimplicialMesh
{

    public:

        SimplicialMesh( int dim, int outerdim );
        virtual ~SimplicialMesh();

        SimplicialMesh( 
            int dim, int outerdim,
            const Coordinates& coords,
            const std::map< std::pair<int,int>, std::vector<IndexMap> >& sub,
            const std::map< std::pair<int,int>, std::vector<std::list<int>> >& super
        );
        
        virtual void check() const;

        virtual void print( std::ostream& out ) const;

        /* Element Access */
        
        int getinnerdimension() const;
        int getouterdimension() const;
        Coordinates& getcoordinates();
        const Coordinates& getcoordinates() const;
        
        int countdimensionscounted() const;
        int countsimplices(int) const;
        const IndexMap getsubsimplices( int, int, int ) const;
        const std::list<int> getsupersimplices( int, int, int ) const;

        /* General management */
        
        bool hassimplexlist(int) const;
        bool hassubsimplexlist(int,int) const;
        bool hassupersimplexlist(int,int) const;
        
        void buildsimplexlist(int);
        void buildsubsimplexlist(int,int);
        void buildsupersimplexlist(int,int);
        
        const std::map< std::pair<int,int>, std::vector<IndexMap> >& getsub() const;
        const std::map< std::pair<int,int>, std::vector<std::list<int>> >& getsuper() const;
        
        void completesimplexlists();
        void completesubsimplexlists();
        void completesupersimplexlists();
        void complete();
        
        /* obtain geometric information */
        
        DenseMatrix getLinearPart( int, int ) const;
        FloatVector getShiftPart( int, int ) const;
        Float getVolume( int, int ) const;
        Float getDiameter( int, int ) const;
        FloatVector getMidpoint( int, int ) const;
        
        
    private:

        int innerdimension;
        int outerdimension;

        Coordinates coordinates;

        /* subcell lists */
        std::map< std::pair<int,int>, std::vector<IndexMap> > subsimplex_list;
        
        /* supercell lists */
        std::map< std::pair<int,int>, std::vector<std::list<int>> > supersimplex_list;
        
};



static inline int countsubsimplices( int n, int k )
{
    assert( 0 <= k && k <= n );
    return binomial<int>( n+1, k+1 );
}




inline std::ostream& operator<<( std::ostream& os, const SimplicialMesh& sm )
{
    sm.print( os );
    return os;
}


#endif