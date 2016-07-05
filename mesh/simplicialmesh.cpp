
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
#include "../operators/floatvector.hpp"
#include "coordinates.hpp"
#include "simplicialmesh.hpp"



SimplicialMesh::SimplicialMesh( int dim, int outerdim )
:
    innerdimension(dim),
    outerdimension(outerdim),
    coordinates(outerdim,0),
    // simplex_list_active(dim+1,false),
    // simplex_list_count(dim+1,0),
    subsimplex_list(),
    supersimplex_list()
{
    // 0 and n dimensional simplices are active, though empty 
    // simplex_list_active[0] = simplex_list_active[dim] = true; // FIXME: New code plan 
    // simplex_list_count[0] = simplex_list_count[dim] = 0;
    
    // there is a n->0 subsimplex list 
    auto key_n_0 = std::make_pair(dim,0);
    subsimplex_list.insert( std::make_pair( key_n_0, std::vector<int>() ) );
}


SimplicialMesh::~SimplicialMesh()
{
    
}


void SimplicialMesh::check() const
{
    // FIXME: add code 
}

void SimplicialMesh::print( std::ostream& out ) const
{
    // FIXME: add code 
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
    // return simplex_list_count.at(d); FIXME: New method
}


const IndexMap SimplicialMesh::getsubsimplices( int dimfrom, int cell, int dimto ) const
{
    assert( 0 <= dimfrom && dimfrom <= getinnerdimension() );
    assert( ! hassimplexlist(dimfrom) );
    assert( 0 <= dimto && dimto < dimfrom );
    assert( ! hassimplexlist(dimto) );
    assert( 0 <= cell && cell < countsimplices(dimfrom) );
    assert( ! hassubsimplexlist(cell,dimto) );
    
    std::vector<int> tempvec = subsimplex_list.find( std::make_pair(cell,dimto) )->second;
    int tempn = countsubsimplices(dimfrom,dimto);
    
    IndexMap ret = IndexMap( 
        IndexRange(0, tempn ), 
        Nat0,
        std::vector<int>( tempvec.begin() + tempn * cell, tempvec.end() + tempn * (cell+1) )
    );
    
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
    // assert( 0 <= d && d <= getinnerdimension() );
    // return simplex_list_active.at( d );
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
	const int counteach = sigmas.size();
	
	/* generate all subsimplices with duplicates */
	std::vector<std::vector<int>> list_full( 
		countsimplices(innerdimension) * counteach,
		std::vector<int>(dim+1) );
	for( int S = 0; S < countsimplices(innerdimension); S++ )
	for( int s = 0; s < counteach; s++ )
	{
		// list_full[ S * counteach + s ]
		// = ( sigmas[ s ] * subsimplex_list[ std::make_pair(innerdimension,0) ].at( S ) ); 
		// FIXME: Check types 
	}
	
	/* eleminate duplicates */
	std::sort( list_full.begin(), list_full.end() );
	std::unique( list_full.begin(), list_full.end() );
	
	/* create new list */
	std::vector<IndexMap>& temp
	=
	subsimplex_list.insert( std::make_pair( std::make_pair(dim,0), std::vector<IndexMap>( list_full.size(), IndexMap( IndexRange(0,dim), IndexRange(0,countsimplices(0)-1) ) ) ) ).first->second;
	
	/* fill up that list */
	for( int t = 0; t < list_full.size(); t++ )
		temp[t] = IndexMap( IndexRange(0,dim), IndexRange(0,countsimplices(0)-1), list_full[t] );
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


    
    
    