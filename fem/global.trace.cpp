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





SparseMatrix FEECBrokenTraceMatrix( const Mesh& mesh, int n, int k, int r, bool is_signed )
{
    
    Assert( 0 <= n );
    Assert( n <= mesh.getinnerdimension() );
    Assert( 0 <= r );
    Assert( 0 <= k and k <= n );

    const std::vector<MultiIndex>& multis = generateMultiIndices( IndexRange( 0, n ), r );
    
    const int dim_polynomials = multis.size();
    // TODO: catch the border cases r <= 0 or k >= n
    /*
    const std::vector<MultiIndex>& multis_src  = generateMultiIndices( IndexRange( 0, n ), r );
    const std::vector<MultiIndex>& multis_dest = generateMultiIndices( IndexRange( 0, n ), r-1 );
    
    const std::vector<IndexMap>& sigmas_cell = generateSigmas( IndexRange( 0, k ), IndexRange( 0, n   ) );
    const std::vector<IndexMap>& sigmas_face = generateSigmas( IndexRange( 0, k ), IndexRange( 0, n-1 ) );
    
    const int dim_sigma_cells = sigmas_cell.size();
    const int dim_sigma_faces = sigmas_face.size();
    
    const int num_cells = mesh.count_simplices(n  );
    const int num_faces = mesh.count_simplices(n-1);

    std::vector< std::vector<IndexMap> > subsimplex_inclusions( n+1 );
    for( auto d : IndexRange(0,n) )
    subsimplex_inclusions[d] = generateSigmas( IndexRange(0,d), IndexRange(0,n) );


    const int dim_in      = num_cells * dim_polynomials * dim_sigma_cells;
    const int dim_out     = num_faces * dim_polynomials * dim_sigma_faces;
    const int num_entries = num_cells * (n+1) * dim_polynomials * dim_sigma_faces;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s <  num_cells; s++ )
    for( int f = 0; f <=         n; f++ )
    {

        for( int mi = 0; mi <      multis.size(); mi++ )
        for( int fs = 0; fs < sigmas_face.size(); fs++ )
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
    
    }

    
    return ret;
    */
    
}



#endif
