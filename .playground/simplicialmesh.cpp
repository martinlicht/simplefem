
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



SimplicialMesh::SimplicialMesh( int innerdim, int outerdim )
:
    innerdimension(innerdim),
    outerdimension(outerdim),
    coordinates(outerdim,0),
    subsimplex_list(),
    supersimplex_list()
{
    // there is a n->0 subsimplex list 
    auto key_n_0 = std::make_pair(innerdim,0);
    subsimplex_list.insert( std::make_pair( key_n_0, std::vector<IndexMap>() ) );
    
    check();
}


SimplicialMesh::SimplicialMesh( int innerdim, int outerdim,
    const Coordinates& coords,
    const std::map< std::pair<int,int>, std::vector<IndexMap> >& sub,
    const std::map< std::pair<int,int>, std::vector<std::list<int>> >& super
)
:
    innerdimension(innerdim),
    outerdimension(outerdim),
    coordinates(coords),
    subsimplex_list(sub),
    supersimplex_list(super)
{
    check();
}


SimplicialMesh::~SimplicialMesh()
{
    
}


void SimplicialMesh::check() const
{
    // FIXME: add code 
    /* checke welche listen vorhanden sein sollen. */
    /* check die alle subsimplex listen. Dabei auch die supersimplex liste pr�fen */
    /* Supersimplex listen checken */
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


int SimplicialMesh::countdimensionscounted() const
{
    int ret = 0;
    for( int d = 0; d <= getinnerdimension(); d++ )
        if( hassimplexlist(d) ) ret++;
    return ret;
}
        

int SimplicialMesh::countsimplices( int d ) const
{
    assert( 0 <= d && d <= getinnerdimension() );
    if( d == 0 ) {
        return coordinates.getnumber();
    } else {
        assert( hassimplexlist(d) );
        return subsimplex_list.find( std::make_pair(d,0) )->second.size();
    }
}


const IndexMap SimplicialMesh::getsubsimplices( int dimfrom, int cell, int dimto ) const
{
    assert( 1 <= dimfrom && dimfrom <= getinnerdimension() );
    assert( hassimplexlist(dimfrom) );
    assert( 0 <= dimto && dimto < dimfrom );
    assert( hassimplexlist(dimto) );
    assert( 0 <= cell && cell < countsimplices(dimfrom) );
    assert( hassubsimplexlist(dimfrom,dimto) );
    
    std::vector<IndexMap> fromto_list = subsimplex_list.find( std::make_pair(dimfrom,dimto) )->second;
    
    IndexMap ret = fromto_list[cell];
    
    return ret;
}


const std::list<int> SimplicialMesh::getsupersimplices( int dimfrom, int cell, int dimto ) const
{
    assert( 1 <= dimto && dimto <= getinnerdimension() );
    assert( hassimplexlist(dimto) );
    assert( 0 <= dimfrom && dimfrom < dimto );
    assert( hassimplexlist(dimfrom) );
    assert( 0 <= cell && cell < countsimplices(dimfrom) );
    assert( hassupersimplexlist(dimfrom,dimto) );
    
    std::vector<std::list<int>> tempvec = supersimplex_list.find( std::make_pair(dimfrom,dimto) )->second;
    
    return tempvec.at( cell );
}



/* General management */

bool SimplicialMesh::hassimplexlist( int d ) const
{
    if( d == 0 )
        return true;
    else 
        return subsimplex_list.end() != subsimplex_list.find(std::make_pair(d,0));
}


bool SimplicialMesh::hassubsimplexlist( int from, int to ) const
{
    assert( 0 <= to && to < from && from <= innerdimension );
    assert( hassimplexlist( from ) && hassimplexlist( to ) );
    auto key = std::make_pair(from,to);
    return subsimplex_list.end() != subsimplex_list.find(key);
}


bool SimplicialMesh::hassupersimplexlist( int from, int to ) const
{
    assert( 0 <= from && from < to && to <= innerdimension );
    assert( hassimplexlist( from ) && hassimplexlist( to ) );
    auto key = std::make_pair(from,to);
    return supersimplex_list.end() != supersimplex_list.find(key);
}



void SimplicialMesh::buildsimplexlist( int dim )
{
    /* set up the basics */
    assert( 0 < dim && dim <= getinnerdimension() );
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
    assert( hassubsimplexlist(from,0) );
    assert( hassubsimplexlist(to,0) );
    assert( ! hassubsimplexlist(from,to) );
    
    const int N = countsubsimplices( from, to );
    auto mynewIMvector = std::vector<IndexMap>( 
        countsimplices(from),
        IndexMap( IndexRange(0,N-1), IndexRange(0,countsimplices(to)-1) )
    );
    
    const std::vector<IndexMap> sigmas = generateSigmas( IndexRange(0,to), IndexRange(0,from) );
    assert( N == sigmas.size() );
    
    for( int S = 0; S < countsimplices(from); S++ )
    for( int s = 0; s < N; s++ )
    {
        int it = 0;
        while( getsubsimplices( to, it, 0 ) != sigmas.at(s) * getsubsimplices( from, S, 0 ) )
            it++;
        assert( it < countsimplices(to) );
        mynewIMvector.at(S)[s] = it;
    }
    
    subsimplex_list.insert( std::make_pair( std::make_pair(from,to), mynewIMvector ) );
    
}


void SimplicialMesh::buildsupersimplexlist( int from, int to )
{
    assert( 0 <= from && from < to && to <= innerdimension );
    assert( hassubsimplexlist(from,0) );
    assert( hassubsimplexlist(to,0) );
    assert( hassubsimplexlist(to,from) );
    assert( ! hassupersimplexlist(from,to) );
    
    std::vector<std::list<int>> mynewlist( countsimplices(from), std::list<int>() );
    
    for( int S = 0; S < countsimplices(to); S++ ) {

        const IndexMap& sub_of_S = getsubsimplices( to, S, from );

        for( int s = 0; s < countsubsimplices(to,from); s++ ) {
            mynewlist[ sub_of_S[s] ].push_back( S );
        }

    }
    
    supersimplex_list.insert( std::make_pair( std::make_pair(from,to), mynewlist ) );
    
}




const std::map< std::pair<int,int>, std::vector<IndexMap> >& SimplicialMesh::getsub() const
{
    return subsimplex_list;
}


const std::map< std::pair<int,int>, std::vector<std::list<int>> >& SimplicialMesh::getsuper() const
{
    return supersimplex_list;
}
        
    
    

void SimplicialMesh::completesimplexlists()
{
    for( int d = 1; d < innerdimension; d++ )
        if( !hassimplexlist( d ) ) buildsimplexlist( d );
}

void SimplicialMesh::completesubsimplexlists()
{
    completesimplexlists();
    for( int D = 1; D < innerdimension; D++ )
        for( int d = 1; d < D; d++ )
            if( !hassubsimplexlist( D, d ) ) buildsubsimplexlist( D, d );
}

void SimplicialMesh::completesupersimplexlists()
{
  completesubsimplexlists();
  for( int D = 1; D < innerdimension; D++ )
        for( int d = 1; d < D; d++ )
            if( !hassupersimplexlist( d, D ) ) buildsupersimplexlist( d, D );
}

void SimplicialMesh::complete()
{
    completesimplexlists();
    completesubsimplexlists();
    completesupersimplexlists();
}
        
            
    
    
    
    
    
    
    
    
    
    
DenseMatrix SimplicialMesh::getLinearPart( int dimension, int cell ) const
{
    assert( 0 <= dimension && dimension <= getinnerdimension() );
    assert( hassimplexlist(dimension) );
    IndexMap im = getsubsimplices( dimension, cell, 0 );
    return getcoordinates().getLinearPart( im );        
}

FloatVector SimplicialMesh::getShiftPart( int dimension, int cell ) const
{
    assert( 0 <= dimension && dimension <= getinnerdimension() );
    assert( hassimplexlist(dimension) );
    IndexMap im = getsubsimplices( dimension, cell, 0 );
    return getcoordinates().getShiftPart( im ); 
}

Float SimplicialMesh::getVolume( int dimension, int cell ) const
{
    assert( 0 <= dimension && dimension <= getinnerdimension() );
    assert( hassimplexlist(dimension) );
    IndexMap im = getsubsimplices( dimension, cell, 0 );
    return Determinant( getcoordinates().getLinearPart( im ) ) / factorial(dimension);
}

Float SimplicialMesh::getDiameter( int dimension, int cell ) const
{
    assert( 0 <= dimension && dimension <= getinnerdimension() );
    assert( hassimplexlist(dimension) );
    IndexMap im = getsubsimplices( dimension, cell, 0 );
    Float ret = 0.;
    for( int i = 0; i <= dimension; i++ )
        for( int j = 0; j <= dimension; j++ )
            if( i != j )
                ret = std::max( ret, ( coordinates.getvectorclone(i) - coordinates.getvectorclone(j) ).norm() );
    return ret;
}

FloatVector SimplicialMesh::getMidpoint( int dimension, int cell ) const
{
    assert( 0 <= dimension && dimension <= getinnerdimension() );
    assert( hassimplexlist(dimension) );
    
    IndexMap im = getsubsimplices( dimension, cell, 0 );
    Float scaling = 1. / (1+getouterdimension());
    FloatVector ret( getouterdimension(), 0. );
    
    for( int i = 0; i <= dimension; i++ )
        ret += coordinates.getvectorclone( im[0], scaling );
    
    return ret;
}


        
        
        
        
        
        
        
        
        
        
        
        
        
        