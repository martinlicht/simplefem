
#include <ostream>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../dense/functions.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/polynomialmassmatrix.hpp"
#include "../fem/global.hodgestarmatrix.hpp"





DenseMatrix EuclideanHodgeStar( int n, int k )
{
    assert( n >= 0 );
    assert( k >= 0 && k <= n );
    
    const std::vector<IndexMap> sigmas_in  = generateSigmas( IndexRange( 1, k   ), IndexRange( 1, n ) );
    const std::vector<IndexMap> sigmas_out = generateSigmas( IndexRange( 1, n-k ), IndexRange( 1, n ) );

    assert( sigmas_in.size() == sigmas_out.size() );

    DenseMatrix ret( sigmas_in.size(), sigmas_out.size(), 0. );

    for( int r = 0; r < sigmas_out.size(); r++ )
    {
        const auto& sigma_out = sigmas_out[r];

        int c = 0;

        for( ; c < sigmas_in.size(); c++ ) 
        {
            const auto& sigma_in = sigmas_in[c];

            bool common_index = false;

            for( auto v : sigma_out.getvalues() )
            for( auto w :  sigma_in.getvalues() )
            {
                common_index = common_index or ( v == w );
            }

            if( not common_index ) break;
        }

        assert( c != sigmas_in.size() );

        const auto& sigma_in = sigmas_in[c];

        int signum = 1.; // determine the specific coefficient needed

        for( int i = 1; i <= k;   i++ )
        for( int j = 1; j <= n-k; j++ )
        {
            if( sigma_in[i] > sigma_out[j] ) signum *= -1;
        }

        ret(r,c) = signum;
    }

    return ret;
}


SparseMatrix FEECBrokenHodgeStarMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    
    const int dim_in      = num_simplices * localdim;
    const int dim_out     = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim;
    
    const std::vector<MultiIndex> multis_in  = generateMultiIndices( IndexRange( 0, n ), r );
    const std::vector<IndexMap>   sigmas_in  = generateSigmas( IndexRange( 1, k   ), IndexRange( 0, n ) );
    
    const std::vector<MultiIndex> multis_out = generateMultiIndices( IndexRange( 0, n ), r );
    const std::vector<IndexMap>   sigmas_out = generateSigmas( IndexRange( 1, n-k ), IndexRange( 0, n ) );
    
    // TODO: complete code 
    // assert( multis_dest.size() * sigmas_dest.size() == multis_src.size()  * sigmas_src.size() );
    assert( dim_out == dim_in );
    
    
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    const DenseMatrix polyMM = polynomialmassmatrix( n, r );

    assert( polyMM.is_square() and polyMM.getdimin() == binomial_integer( n+r, n ) );
    
    /*
    Method 1:

    1. Develop a QR decomposition of the transformation matrix to obtain local Euclidean coordinate system;
       In the full-dimensional case, an alternative method is feasible as well (introduce boolean switch).
    2. Write barycentric derivatives in terms of those coordinates (canonicalize the zeroth gradient)
    3. Apply Euclidean Hodge star
    4. Transform back into barycentric coordinate system 

    Method 2:
    1. Write the mass pairing of the n-k forms 
    2. Define the Hodge star implicitly.
    */

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        // obtain geometric information 
        
        Float measure      = mesh.getMeasure( n, s );

        DenseMatrix GPM    = mesh.getGradientProductMatrix( n, s );
          
        assert( measure >= 0. );

        // create local mass matrix 
        
        DenseMatrix formMM = SubdeterminantMatrix( GPM, k );
    
        DenseMatrix fullMM = MatrixTensorProduct( polyMM, formMM ) * measure;

        // create local Hodge product matrix 
        
        DenseMatrix formHM(0); // TODO: wierd ...

        DenseMatrix fullHM = MatrixTensorProduct( polyMM, formHM ) * measure;

        // compute the final local matrix and register its entries 

        DenseMatrix localmatrix = Inverse(fullMM) * fullHM;

        for( int i = 0; i < localdim; i++ )
        for( int j = 0; j < localdim; j++ )
        {
            int index_of_entry = s * localdim * localdim + i * localdim + j;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * localdim + i;
            entry.column = s * localdim + j;
            entry.value  = localmatrix( i, j );
            
            ret.setentry( index_of_entry, entry );
        }
        
    }
    
    return ret;
}


