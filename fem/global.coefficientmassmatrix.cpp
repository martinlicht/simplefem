
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/cholesky.hpp"
#include "../dense/functions.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/matrixtensorproduct.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/local.polynomialmassmatrix.hpp"

#include "../fem/global.massmatrix.hpp"



SparseMatrix FEECBrokenCoefficientMassMatrix( const Mesh& mesh, int n, int k, int r,
                                              int s, std::function<DenseMatrix(FloatVector)>& generator 
) {
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    assert( s >= 0 );
    
    // Auxiliary calculations and preparations
    
    // generate the polynomial mass matrices wrt to same base weight 
    // create the coefficients of the lagrange polynomials for each integration point 
    // sum up so that you get the polynomial mass matrices for each integration point 

    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    
    const int internaldim = binomial_integer( n+s, n )

    const int dim_in      = num_simplices * localdim;
    const int dim_out     = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim * internaldim;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    auto internal_mi = generateMultiIndices( IndexRange(0,n), s );

    std::vector<DenseMatrix> coefficientmassmatrices( internaldim );
    
    for( int i = 0; i < coefficientmassmatrices.size(); i++ )
        coefficientmassmatrices[i]
    
    DenseMatrix polyMM = polynomialmassmatrix( n, r ); // TODO: polymassmatrix with extra base weight 

    assert( polyMM.issquare() and polyMM.getdimin() == binomial_integer( n+r, n ) );
    
//     LOG << polyMM << nl;
        
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        /* get the values of form matrix at the integration points */
        /* get the polynomial coefficients of the corresponding Lagrange polynomial */


        Float measure       = mesh.getMeasure( n, s );

        DenseMatrix GM      = mesh.getGradientMatrix( n, s );
        
        DenseMatrix Product = GM // TODO: complete this
        assert(false);
        DenseMatrix formMM  = SubdeterminantMatrix( Product, k );
    
        DenseMatrix fullMM  = MatrixTensorProduct( polyMM, formMM ) * measure;

        assert( measure >= 0. );

        // LOG << measure << nl;
        
        // LOG << formMM << nl;
        
        // LOG << fullMM << nl;
        
        for( int i = 0; i < localdim; i++ )
        for( int j = 0; j < localdim; j++ )
        {
            int index_of_entry = s * localdim * localdim + i * localdim + j;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * localdim + i;
            entry.column = s * localdim + j;
            entry.value  = fullMM( i, j );
            
            ret.setentry( index_of_entry, entry );
        }
        
        
        
    }
    
    return ret;
}




