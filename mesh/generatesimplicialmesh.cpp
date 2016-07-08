
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






SimplicialMesh UnitCubeTriangulation( int innerdim, int outerdim )
{

    // TODO: Kuhn triangulatin of unit cube 
    // * Lade die Liste der Permutationen
    // * Rename this method to KuhnTriangulation
    // * make it a static method
    // Lade die 2^n verschiedenen vertices 
    // Erstelle die n simplices mittels iteration ueber die Permutationen
    
    assert( innerdim <= outerdim     );
    
    Coordinates coords( outerdim, 1 << innerdim );
    for( int v = 0; v < ( 1 << innerdim ); v++ ) {
        for( int p = 0; p < innerdim; p++ )
            coords.getvector(v).setentry( p, ( v >> p ) % 2 );
        for( int p = innerdim; p < outerdim; p++ )
            coords.getvector(v).setentry( p, 0. );
    }
        
    IndexRange ir( 0, innerdim-1 );
    std::vector<IndexMap> perms = generatePermutations( ir );
    
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

        for( int p = 1; p < innerdim+1; p++ ) {
            index += 1 << (*pi)[p];
            newsimplex[p] = index;
        }
        
        newsimplices[picounter] = IndexMap( IndexRange(0,innerdim),
                                            IndexRange(0,coords.getnumber()-1),
                                            newsimplex );
    }
    
    // FIXME: Koordinaten und Simplizes verarbeiten 
    
}


void generateMeshFromStream( std::istream& in )
{
    // TODO: add code 
}

