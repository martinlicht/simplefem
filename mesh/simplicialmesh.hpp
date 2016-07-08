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
#include "coordinates.hpp"




class SimplicialMesh
{

public:

    SimplicialMesh( int dim, int outerdim );
    virtual ~SimplicialMesh();

    SimplicialMesh( 
        int dim, int outerdim,
        const Coordinates& coords,
        std::map< std::pair<int,int>, std::vector<IndexMap> >& sub,
        std::map< std::pair<int,int>, std::vector<std::list<int>> >& super 
    );
    
    virtual void check() const;

    virtual void print( std::ostream& out ) const;

    /* Element Access */
    
    int getinnerdimension() const;
    int getouterdimension() const;
    Coordinates& getcoordinates();
    const Coordinates& getcoordinates() const;
    
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
    
    /* Construction of a concrete mesh */
    
    void addfromstream( std::istream& in );
    
    void addfrom( const SimplicialMesh& );
    
    void addunitcube( const FloatVector&, Float );
        
private:

    int innerdimension;
    int outerdimension;

    Coordinates coordinates;

    /* subcell lists */
    std::map< std::pair<int,int>, std::vector<IndexMap> > subsimplex_list;
    
    /* supercell lists */
    std::map< std::pair<int,int>, std::vector<std::list<int>> > supersimplex_list;
    
};



inline int countsubsimplices( int n, int k )
{
    assert( 0 <= k && k <= n );
    return binomial( n+1, k+1 );
}






#endif