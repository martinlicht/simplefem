
#include <vector>
#include <iostream>
#include <algorithm>
// #include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/densematrix.hpp"
#include "../operators/linearoperator.hpp"

#include "element.matrix.hpp"









int polysigmaindex2fullindex( int n, int r, int k, int poly_i, int sigma_i )
{
  int poly_n = binomial( n + r, r );
  int sigma_n = binomial( n + 1, k );
  assert( 0 <= poly_i && poly_i < poly_n );
  assert( 0 <= sigma_i && sigma_i < sigma_n );
  return poly_i * poly_n + sigma_i;
}

int fullindex2polysigmaindex( int n, int r, int k, int full_i, int& poly_i, int& sigma_i )
{
  int poly_n = binomial( n + r, r );
  int sigma_n = binomial( n + 1, k );
  assert( 0 <= full_i && full_i < poly_n * sigma_n );
  poly_i = full_i / poly_n;
  sigma_i = full_i % poly_n;
}









DenseMatrix calculateElementDiffMatrix(
                int n, // innerdimension
                int r, // polynomial degree
                int k // form degree 
                )
{
    assert( 0 <= n );
    assert( 0 <= r );
    assert( 0 <= k ); //FIXME: Decide whether to allow larger index range.
    
    IndexRange srcrange = IndexRange( 0, k );
    IndexRange tarrange = IndexRange( 0, k+1 );
    IndexRange coordrange = IndexRange( 0, n );
    
    std::vector<IndexMap> sigmas_src = generateSigmas( srcrange, coordrange );
    std::vector<IndexMap> sigmas_tar = generateSigmas( tarrange, coordrange );
    
    std::vector<MultiIndex> multis_src = generateMultiIndices( coordrange, r );
    std::vector<MultiIndex> multis_tar = generateMultiIndices( coordrange, r );
    
    DenseMatrix ret( sigmas_tar.size() * multis_tar.size(), 
                     sigmas_src.size() * multis_src.size(), 
                     0 );
    
    for( int sigma_src_i = 0; sigma_src_i < sigmas_src.size(); sigma_src_i++ )
    for( int poly_src_i = 0; poly_src_i < multis_src.size(); poly_src_i++ )
    for( const auto c : coordrange )
    {
        
        /* obtain target multiindex and factor  */
        const MultiIndex& poly_src = multis_src[poly_src_i];
        
        if( poly_src[c] == 0 )
          continue;
          
        MultiIndex poly_tar = poly_src - c;
        
        int factor = poly_src[c];
        
        /* obtain target sigma and the signum */
        
        const IndexMap& sigma_src = sigmas_src[sigma_src_i];
        
        if( sigma_src.rangecontains( c ) )
          continue;
        
        IndexMap sigma_tar( tarrange, coordrange );
        for( const auto& sigma_test : sigmas_tar )
          if( sigma_test.rangecontains(c) ){
            bool flag = true;
            for( const auto& i : srcrange )
              flag &= sigma_test.rangecontains( sigma_src[i] );
            if(flag)
              sigma_tar = sigma_test;
          }
        
        int signum = integerpower( -1, sigma_tar.preimageof(c) );
        
        /* find positions */
        int poly_tar_i = find( multis_tar.begin(), multis_tar.end(), poly_tar ) - multis_tar.begin();
        int sigma_tar_i = find( sigmas_tar.begin(), sigmas_tar.end(), sigma_tar ) - sigmas_tar.begin();
        assert( poly_tar_i < multis_tar.size() );
        assert( sigma_tar_i < sigmas_tar.size() );
        
        /* set matrix entry */
        int row_index = polysigmaindex2fullindex( n, r  , k  , poly_src_i, sigma_src_i ); 
        int col_index = polysigmaindex2fullindex( n, r-1, k+1, poly_tar_i, sigma_tar_i ); 
        ret( row_index, col_index ) = signum * factor;
    }
    
    return ret;
    
}



