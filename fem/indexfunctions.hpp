#ifndef INCLUDEGUARD_FEM_INDEXFUNCTIONS
#define INCLUDEGUARD_FEM_INDEXFUNCTIONS


#include <iostream>
#include <utility>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"




//////////////////////////////////////////////////////
//                                                  //
//  List of Sullivan Indices                        //
//                                                  //
//  alpha is a multiindex                           //
//  sigma 1:k -> 0:n                                //
//  min[alpha] notin [sigma]                        //
//  [alpha] u [sigma] = [0..n]                      //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////




inline std::vector< std::pair<MultiIndex,IndexMap> > ListOfSullivanIndices( int n, int k, int r )
{
    
    // check whether the parameters are right 
    
    assert( r >= 1 );
    assert( n >= 0 );
    assert( k >= 0 );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    
    if( k > n ) 
        return std::vector< std::pair<MultiIndex,IndexMap> >();
    
    // Auxiliary calculations and preparations
    
//     const IndexRange N( 0, n );
    const std::vector<int> N = [&n]()->auto{ std::vector<int> ret(n+1); 
                                           for( int i = 0; i <= n; i++ ) ret[i] = i;
                                           return ret; }();
    
    const std::vector<MultiIndex> alphas = generateMultiIndices( IndexRange( 0, n ), r );
    const std::vector<IndexMap>   sigmas = generateSigmas( IndexRange( 1, k ), IndexRange( 0, n ) );
    
    std::vector< std::pair<MultiIndex,IndexMap> > ret;
    
    //  [ size of set taken from Acta paper ]
    int computed_length = binomial_integer( r-1, n-k ) * binomial_integer( r+k, k );
    
    for( const MultiIndex& alpha : alphas )
    for( const IndexMap&   sigma : sigmas )
    {
        
        // First, check that every p in 0..n is contained in the ranges of alpha and/or sigma
        bool b1 = std::all_of( N.begin(), N.end(), 
                               [ &alpha, &sigma ]( int p ) -> bool { 
                                   return sigma.rangecontains(p) or alpha[p] > 0;
                                   }
                             );
        
        // Second, check that min[alpha] not in [sigma]
        bool b2 = not sigma.rangecontains( alpha.min() );
        
        // if both criteria are satisfied, then save that one
        if( b1 and b2 )
            ret.push_back( std::pair<MultiIndex,IndexMap>( alpha, sigma ) );
        
    }
    
    assert( ret.size() == computed_length );
    
    return ret;
    
}










#endif
