
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
    assert( coords.getdimension() == outerdim );
    assert( coords.getnumber() == integerpower(2,innerdim) );
    cout << "############" << coords.getnumber() << std::endl;
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



SimplicialMesh UnitCubeTriangulation2D( int xstep, int ystep )
{
    assert( 1 <= xstep && 1 <= ystep );
    
    Coordinates coords( 3, (xstep+1)*(ystep+1) );
    for( int v = 0; v < coords.getnumber(); v++ ) {
      
        Float xcoord = ( v / (ystep+1) ) / (Float)xstep;
        Float ycoord = ( v % (ystep+1) ) / (Float)ystep;
        Float zcoord = xcoord * ycoord;
        
        coords.setdata( v, 0, xcoord );
        coords.setdata( v, 1, ycoord );
        coords.setdata( v, 2, zcoord );
    }
    
    std::vector<IndexMap> newsimplices( 2 * xstep * ystep, 
                                    IndexMap( IndexRange(0,2), 
                                              IndexRange(0,coords.getnumber()-1)
                                            )
                                    );

    for( int x = 0; x < xstep; x++ )
    for( int y = 0; y < ystep; y++ )
    {
        int base = 2 * ( x * ystep + y );
        assert( base + 1 < 2 * xstep * ystep );
        newsimplices[base+0][0] = ( x + 0 ) * (ystep+1) + y + 0;
        newsimplices[base+0][1] = ( x + 1 ) * (ystep+1) + y + 0;
        newsimplices[base+0][2] = ( x + 1 ) * (ystep+1) + y + 1;
        newsimplices[base+1][0] = ( x + 0 ) * (ystep+1) + y + 0;
        newsimplices[base+1][1] = ( x + 0 ) * (ystep+1) + y + 1;
        newsimplices[base+1][2] = ( x + 1 ) * (ystep+1) + y + 1;
    }
        
    std::map< std::pair<int,int>, std::vector<IndexMap> > new_subsimplexlist;
    new_subsimplexlist.insert( std::make_pair( std::make_pair(2,0), newsimplices ) );
    
    std::map< std::pair<int,int>, std::vector<std::list<int>> > dummy_supersimplexlist;
    
    return SimplicialMesh( 2, 3, coords, new_subsimplexlist, dummy_supersimplexlist );
    
}



SimplicialMesh UnitCubeTriangulation3D( int xstep, int ystep, int zstep )
{
    assert( 1 <= xstep && 1 <= ystep && 1 <= zstep );
    
      Coordinates coords( 3, (xstep+1)*(ystep+1)*(zstep+1) );
    for( int x = 0; x <= xstep; x++ )
    for( int y = 0; y <= ystep; y++ )
    for( int z = 0; z <= zstep; z++ )
    {
        Float xcoord = x / (Float)xstep;
        Float ycoord = y / (Float)ystep;
        Float zcoord = z / (Float)zstep;
        
        int v = x * (ystep+1)*(zstep+1) + y * (zstep+1) + z;
        assert( 0 <= v && v < coords.getnumber() );
        
        coords.setdata( v, 0, xcoord );
        coords.setdata( v, 1, ycoord );
        coords.setdata( v, 2, zcoord );
    }
    
    std::vector<IndexMap> newsimplices( 6 * xstep * ystep * zstep, 
                                        IndexMap( IndexRange(0,3), 
                                                  IndexRange(0,coords.getnumber()-1)
                                                )
                                        );
    
    for( int x = 0; x < xstep; x++ )
    for( int y = 0; y < ystep; y++ )
    for( int z = 0; z < zstep; z++ )
    {
    
        int base = 6 * ( x * ystep * zstep + y * zstep + z );
        assert( base + 5 < 6 * xstep * ystep * zstep );
        // Simplex z y x 
        newsimplices.at(base+0)[0] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 0;
        newsimplices.at(base+0)[1] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 1;
        newsimplices.at(base+0)[2] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 1;
        newsimplices.at(base+0)[3] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 1;
        // Simplex z x y 
        newsimplices.at(base+1)[0] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 0;
        newsimplices.at(base+1)[1] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 1;
        newsimplices.at(base+1)[2] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 1;
        newsimplices.at(base+1)[3] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 1;
        // Simplex x z y 
        newsimplices.at(base+2)[0] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 0;
        newsimplices.at(base+2)[1] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 0;
        newsimplices.at(base+2)[2] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 1;
        newsimplices.at(base+2)[3] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 1;
        // Simplex x y z 
        newsimplices.at(base+3)[0] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 0;
        newsimplices.at(base+3)[1] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 0;
        newsimplices.at(base+3)[2] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 0;
        newsimplices.at(base+3)[3] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 1;
        // Simplex y x z 
        newsimplices.at(base+4)[0] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 0;
        newsimplices.at(base+4)[1] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 0;
        newsimplices.at(base+4)[2] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 0;
        newsimplices.at(base+4)[3] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 1;
        // Simplex y z x  
        newsimplices.at(base+5)[0] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 0 ) * (zstep+1) + z + 0;
        newsimplices.at(base+5)[1] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 0;
        newsimplices.at(base+5)[2] = ( x + 0 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 1;
        newsimplices.at(base+5)[3] = ( x + 1 ) * (ystep+1) * (zstep+1) + ( y + 1 ) * (zstep+1) + z + 1;
        
    }
        
    std::map< std::pair<int,int>, std::vector<IndexMap> > new_subsimplexlist;
    new_subsimplexlist.insert( std::make_pair( std::make_pair(3,0), newsimplices ) );
    
    std::map< std::pair<int,int>, std::vector<std::list<int>> > dummy_supersimplexlist;
    
    return SimplicialMesh( 3, 3, coords, new_subsimplexlist, dummy_supersimplexlist );
    
}







void generateMeshFromStream( std::istream& )
{
    // TODO: add code 
}

