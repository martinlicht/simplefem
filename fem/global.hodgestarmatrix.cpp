
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


