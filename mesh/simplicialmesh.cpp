
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>


#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "coordinates.hpp"
#include "simplicialmesh.hpp"



SimplicialMesh::SimplicialMesh( int dim, int outerdim )
:
    innerdimension(dim),
    outerdimension(outerdim),
    coordinates(outerdim,0),
    subsimplex_list(),
    supersimplex_list()
{
    // FIXME: Write constructor 
    
    // there is a n->0 subsimplex list 
    auto key_n_0 = std::make_pair(dim,0);
    subsimplex_list.insert( std::make_pair( key_n_0, std::vector<IndexMap>() ) );
}


SimplicialMesh::~SimplicialMesh()
{
    
}


void SimplicialMesh::check() const
{
    // FIXME: add code 
}

void SimplicialMesh::print( std::ostream& os ) const
{
    os << "Printe Mesh!" << std::endl;
}


/* Element Access */

int SimplicialMesh::getinnerdimension() const
{
    return innerdimension;
}

int SimplicialMesh::getouterdimension() const
{
    return outerdimension;
}

Coordinates& SimplicialMesh::getcoordinates()
{
    return coordinates;
}


const Coordinates& SimplicialMesh::getcoordinates() const
{
    return coordinates;
}



int SimplicialMesh::countsimplices( int d ) const
{
    assert( 0 <= d && d <= getinnerdimension() );
    assert( hassimplexlist(d) );
    return subsimplex_list.find( std::make_pair(d,0) )->second.size();
}


const IndexMap SimplicialMesh::getsubsimplices( int dimfrom, int cell, int dimto ) const
{
    assert( 0 <= dimfrom && dimfrom <= getinnerdimension() );
    assert( ! hassimplexlist(dimfrom) );
    assert( 0 <= dimto && dimto < dimfrom );
    assert( ! hassimplexlist(dimto) );
    assert( 0 <= cell && cell < countsimplices(dimfrom) );
    assert( ! hassubsimplexlist(cell,dimto) );
    
    std::vector<IndexMap> fromto_list = subsimplex_list.find( std::make_pair(cell,dimto) )->second;
    
    IndexMap ret = fromto_list[cell];
    
    return ret;
}


const IndexMap SimplicialMesh::getsupersimplices( int dimfrom, int cell, int dimto ) const
{
    assert( 0 <= dimto && dimto <= getinnerdimension() );
    assert( ! hassimplexlist(dimto) );
    assert( 0 <= dimfrom && dimfrom < dimto );
    assert( ! hassimplexlist(dimfrom) );
    assert( 0 <= cell && cell < countsimplices(dimfrom) );
    assert( ! hassupersimplexlist(dimfrom,dimto) );
    
    std::vector<IndexMap> tempvec = supersimplex_list.find( std::make_pair(dimfrom,dimto) )->second;
    
    return tempvec.at( cell );
}



/* General management */

bool SimplicialMesh::hassimplexlist( int d ) const
{
    return hassubsimplexlist( d, 0 );
}


bool SimplicialMesh::hassubsimplexlist( int from, int to ) const
{
    assert( 0 <= to && to < from && from <= innerdimension );
    auto key = std::make_pair(from,to);
    return subsimplex_list.end() != subsimplex_list.find(key);
}


bool SimplicialMesh::hassupersimplexlist( int from, int to ) const
{
    assert( 0 <= from && from < to && to <= innerdimension );
    auto key = std::make_pair(from,to);
    return supersimplex_list.end() != supersimplex_list.find(key);
}



void SimplicialMesh::buildsimplexlist( int dim )
{
    /* set up the basics */
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( ! hassimplexlist(dim) );
    std::vector<IndexMap> sigmas
    = generateSigmas( IndexRange(0,dim), IndexRange(0,innerdimension) );
    const int countsigmas = sigmas.size();

    /* generate all subsimplices with duplicates */
    std::vector<IndexMap> list_full( 
            countsimplices(innerdimension) * countsigmas,
            IndexMap( IndexRange(0,dim), IndexRange(0,countsimplices(0)) ) );

    for( int S = 0; S < countsimplices(innerdimension); S++ )
    for( int s = 0; s < countsigmas; s++ )
    {
        list_full[ S * countsigmas + s ]
        = sigmas[ s ] * subsimplex_list[ std::make_pair(innerdimension,0) ].at( S ); 
    }

    /* eleminate duplicates */
    std::sort( list_full.begin(), list_full.end() );
    std::unique( list_full.begin(), list_full.end() );
    
    subsimplex_list.insert( std::make_pair( std::make_pair(dim,0), list_full ) );
}
    

void SimplicialMesh::buildsubsimplexlist( int from, int to )
{
    assert( 0 <= to && to < from && from <= innerdimension );
    assert( ! hassubsimplexlist(from,to) );
    // FIXME: add code 
}


void SimplicialMesh::buildsupersimplexlist( int from, int to )
{
    assert( 0 <= from && from < to && to <= innerdimension );
    assert( ! hassupersimplexlist(from,to) );
    // FIXME: add code 
}


    
    
/* Construction of a concrete mesh */

void SimplicialMesh::addfromstream( std::istream& in )
{
    // FIXME: add code 
}

void SimplicialMesh::addfrom( const SimplicialMesh& )
{
    // FIXME: add code 
}

void SimplicialMesh::addunitcube( const FloatVector&, Float )
{
    // FIXME: add code 
}


    
    
    