
#include "generatesimplicialmesh.hpp"

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>


#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "coordinates.hpp"
#include "simplicialmesh.hpp"

using namespace std;




SimplicialMesh UnitCubeTriangulation( int innerdim, int outerdim )
{

    // NOTE: Kuhn triangulatin of unit cube 
    // * Lade die Liste der Permutationen
    // * Rename this method to KuhnTriangulation
    // * make it a static method
    // Lade die 2^n verschiedenen vertices 
    // Erstelle die n simplices mittels iteration ueber die Permutationen
    
    assert( innerdim <= outerdim );
    
	Coordinates coords( outerdim, integerpower(2,innerdim) );
    for( int v = 0; v < coords.getnumber(); v++ ) {
        for( int p = 0; p < innerdim; p++ ) {
			coords.setdata( v, p, getbit(v,p) );
		}
        for( int p = innerdim; p < outerdim; p++ )
            coords.setdata( v, p, 0. );
    }
    
	IndexRange ir( 1, innerdim );
    std::vector<IndexMap> perms = generatePermutations( ir ); // generateIndexMaps( ir, ir );
	
	std::vector<IndexMap> newsimplices( perms.size(), 
                                        IndexMap( IndexRange(0,innerdim), 
                                                  IndexRange(0,coords.getnumber()-1)
                                                )
                                        );
    
	int picounter=0;
    for( auto pi = perms.begin(); pi != perms.end(); pi++, picounter++ )
    {
        int index = 0;
		
		std::vector<int> newsimplex(innerdim+1);
        newsimplex[0] = index;

        for( int p = 1; p <= innerdim; p++ ) {
            index += integerpower( 2, (*pi)[p]-1 );
            newsimplex[p] = index;
			assert( index <= coords.getnumber()-1 );
        }
        
		newsimplices[picounter] = IndexMap( IndexRange(0,innerdim),
                                            IndexRange(0,coords.getnumber()-1),
                                            newsimplex );
    }
    
    std::map< std::pair<int,int>, std::vector<IndexMap> > thatthing;
    thatthing.insert( std::make_pair( std::make_pair(innerdim,0), newsimplices ) );
    
    std::map< std::pair<int,int>, std::vector<std::list<int>> > thatotherthing;
    
    return SimplicialMesh( innerdim, outerdim, coords, thatthing, thatotherthing );
    
}


void generateMeshFromStream( std::istream& in )
{
    // TODO: add code 
}

