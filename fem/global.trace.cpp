#ifndef INCLUDEGUARD_FEM_TRACEMATRIX
#define INCLUDEGUARD_FEM_TRACEMATRIX


#include <vector>
#include <iostream>
// #include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"





SparseMatrix BrokenTraceMatrix( const Mesh& mesh, int n, int k, int r, bool is_signed )
{
    
    // TODO: catch the border cases r <= 0 or k >= n
    
    const std::vector<MultiIndex>& multis_src  = generateMultiIndices( IndexRange( 0, n ), r );
    const std::vector<MultiIndex>& multis_dest = generateMultiIndices( IndexRange( 0, n ), r-1 );
    
    const std::vector<IndexMap>&   sigmas_src  = generateSigmas( IndexRange( 0, k ), IndexRange( 0, n ) );
    const std::vector<IndexMap>&   sigmas_dest = generateSigmas( IndexRange( 0, k+1 ), IndexRange( 0, n ) );
    
    DenseMatrix ret( multis_dest.size() * sigmas_dest.size(), multis_src.size() * sigmas_scr.size(), 0. )
    
    for( int i = 0; i < multis_src.size(); i++ )
    for( int j = 0; j < sigmas_src.size(); j++ )
    for( int c = 0; c <= n; c++ )
    {
        
        const MultiIndex& multi_src = multis_src[i];
        const MultiIndex& sigma_src = sigmas_src[j];
        
        const int index_col = i * multis_src.size() + j;
        
        if( multis_src[i][c] == 0 ) continue;
        if( sigmas_src[j].rangecontains( c ) ) continue;
        
        MultiIndex multi_dest = multi_src - c;
        MultiIndex sigma_dest = sigma_src + c;
        
        int k = std::find( multis_dest.begin(), multis_dest.end(), multi_dest ) - multis_dest.begin();
        int l = std::find( multis_dest.begin(), multis_dest.end(), sigma_dest ) - sigmas_dest.begin();
        
        const int index_row = k * multis_src.size() + l;
        
        ret( index_row, index_col ) = multi_src[c] * sign( c, sigma_src );
        
    }
    
    return ret;
    
}



#endif
