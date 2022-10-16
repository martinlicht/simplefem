
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
#include "../fem/utilities.hpp"

#include "../fem/local.polynomialmassmatrix.hpp"

#include "../fem/global.coefficientmassmatrix.hpp"



SparseMatrix FEECBrokenCoefficientMassMatrix( const Mesh& mesh, int n, int k, int r,
                                              int w, const std::function<DenseMatrix(const FloatVector&)>& generator 
) {
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    assert( w >= 0 );
    
    // Dimensions of the output matrix and number of entries 
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    
    const int dim_in      = num_simplices * localdim;
    const int dim_out     = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    
    // assemble algebraic auxiliary material
    // - lagrange points in barycentric coordinates 
    // - coefficients of Lagrange polynomials
    // mass matrices 

    const auto lpbcs = InterpolationPointsInBarycentricCoordinates( n, w );

    const auto lpcoeff = Inverse( PointValuesOfMonomials( w, lpbcs ) );
    
    const auto polymassmatrix_per_point = polynomialmassmatrices_per_lagrangepoint( n, r, w );
    
    
    // loop over the simplices and compute the mass matrices

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        // assemble some data for the element 
        // - measure 
        // - barycentric coordinates 
        // - lagrange points 
        
        Float measure       = mesh.getMeasure( n, s );
        assert( measure >= 0. );

        DenseMatrix GM    = mesh.getGradientMatrix( n, s );
        DenseMatrix extGM = SubdeterminantMatrix( GM, k );

        auto vertex_coordinates = mesh.getVertexCoordinateMatrix( n, s );
        auto lpeucl             = vertex_coordinates * lpbcs;

        // compute the mass matrix contribution 
        // for each lagrange point 
        
        DenseMatrix full_element_matrix( localdim, localdim, 0. );

        for( int p = 0; p < polymassmatrix_per_point.size(); p++ )
        {
            DenseMatrix matrix_at_point = generator( lpeucl.getcolumn(p) );

            auto polyMM = polymassmatrix_per_point[p];
            
            auto formMM = Transpose(extGM) * matrix_at_point * extGM;

            // DenseMatrix GPM = SubdeterminantMatrix( mesh.getGradientProductMatrix( n, s ), k );
            // assert( ( GPM - formMM ).issmall() );

            if( w == 0 ) assert( ( polyMM - polynomialmassmatrix(n,r) ).issmall() );

            auto fullMM = measure * MatrixTensorProduct( polyMM, formMM );

            full_element_matrix = full_element_matrix + fullMM;
        }
        
        // DONE ... now list everything.

        for( int row = 0; row < localdim; row++ )
        for( int col = 0; col < localdim; col++ )
        {
            int index_of_entry = s * localdim * localdim + row * localdim + col;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * localdim + row;
            entry.column = s * localdim + col;
            entry.value  = full_element_matrix( row, col );
            
            ret.setentry( index_of_entry, entry );
        }

    }

    LOG << "Finished Sparse Matrix entries\n";
    
    return ret;
}




