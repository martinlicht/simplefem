
#include <ostream>
#include <vector>

#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/cholesky.hpp"
#include "../dense/functions.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/matrixtensorproduct.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/local.polynomialmassmatrix.hpp"

#include "../fem/global.massmatrix.hpp"



SparseMatrix FEECBrokenMassMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
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
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    DenseMatrix polyMM = polynomialmassmatrix( n, r );

    assert( polyMM.issquare() and polyMM.getdimin() == binomial_integer( n+r, n ) );
    
//     LOG << polyMM << std::endl;
        
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        Float measure      = mesh.getMeasure( n, s );

        DenseMatrix GPM    = mesh.getGradientProductMatrix( n, s );
            
        DenseMatrix formMM = SubdeterminantMatrix( GPM, k );
    
        DenseMatrix fullMM = MatrixTensorProduct( polyMM, formMM ) * measure;

        assert( measure >= 0. );

        if( k == 0 )
        {
            assert( ( fullMM - polyMM * measure ).issmall() );
        }
        
        // LOG << measure << std::endl;
        
        // LOG << formMM << std::endl;
        
        // LOG << fullMM << std::endl;
        
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





SparseMatrix FEECBrokenMassMatrixRightFactor( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim_in  = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    const int localdim_out = binomial_integer( n+r, n ) * binomial_integer( n  , k );

    const int dim_in      = num_simplices * localdim_in;
    const int dim_out     = num_simplices * localdim_out;
    const int num_entries = num_simplices * localdim_in * localdim_out;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    DenseMatrix polyMM = polynomialmassmatrix( n, r );
    
    DenseMatrix polyMM_right = Transpose(CholeskyDecomposition(polyMM));

    {

        // LOG << polyMM << std::endl;        
        // LOG << ( Transpose(polyMM_right) * polyMM_right ) << std::endl;        
        // LOG << polyMM_right.getdimin() << std::endl;        
        // LOG << polyMM.getdimin() << space << ( Transpose(polyMM_right) * polyMM_right - polyMM ).norm() << std::endl;        
        assert( ( Transpose(polyMM_right) * polyMM_right - polyMM ).issmall() ); 

    }
            
//     LOG << polyMM << std::endl;
        
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        Float measure = mesh.getMeasure( n, s );
            
        DenseMatrix GPM_right    = mesh.getGradientProductMatrixRightFactor( n, s );
            
        DenseMatrix formMM_right = SubdeterminantMatrix( GPM_right, k );
    
        DenseMatrix fullMM_right = MatrixTensorProduct( polyMM_right, formMM_right ) * std::sqrt(measure);

        assert( fullMM_right.getdimin()  == formMM_right.getdimin()  * polyMM_right.getdimin()  );
        assert( fullMM_right.getdimout() == formMM_right.getdimout() * polyMM_right.getdimout() );
        
        {
            
            /* CHECK WHETHER THE ARITHMETICS WORK OUT */
            
            DenseMatrix polyMM_0 = polynomialmassmatrix( n, r );
            DenseMatrix formMM_0 = SubdeterminantMatrix( mesh.getGradientProductMatrix( n, s ), k );
            DenseMatrix fullMM_0 = MatrixTensorProduct( polyMM_0, formMM_0 ) * measure;
            
            assert( ( Transpose(formMM_right) * formMM_right - formMM_0 ).issmall() ); 
            assert( ( Transpose(fullMM_right) * fullMM_right - fullMM_0 ).issmall() ); 
            
        }
        
        
        for( int i = 0; i < localdim_out; i++ )
        for( int j = 0; j < localdim_in ; j++ )
        {
            int index_of_entry = s * localdim_out * localdim_in + i * localdim_in + j;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * localdim_out + i;
            entry.column = s * localdim_in  + j;
            entry.value  = fullMM_right( i, j );
            
            ret.setentry( index_of_entry, entry );
        }
        
        
        
    }
    
    return ret;
}





FloatVector FEECBrokenMassMatrix_cellwisemass( const Mesh& mesh, int n, int k, int r, const FloatVector vec )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
        
    FloatVector ret( num_simplices );
    
    assert( vec.getdimension() == localdim * num_simplices );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        ret[s] = 0.;
        for( int i = 0; i < localdim; i++ )
            ret[s] = ret[s] + vec[ s * localdim + i] * vec[ s * localdim + i];
        
    }
    
    return ret;
}



SparseMatrix FEECBrokenMassMatrix_cellwiseinverse( const Mesh& mesh, int n, int k, int r )
{
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    
    const int dim_in      = num_simplices * localdim;
    const int dim_out     = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    DenseMatrix polyMM = polynomialmassmatrix( n, r );
    
    DenseMatrix polyMMinv = Inverse( polyMM ); 
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        Float measure      = mesh.getMeasure( n, s );

        DenseMatrix GPM    = mesh.getGradientProductMatrix( n, s );
        
        auto aux1   = GPM + DenseMatrix( n+1, n+1, 1. );
        auto aux2   = Inverse( aux1 );
        auto GPMinv = aux2 - DenseMatrix( n+1, n+1, 1. );
            
        DenseMatrix formMMinv = SubdeterminantMatrix( GPMinv, k );
    
        DenseMatrix fullMMinv = MatrixTensorProduct( polyMMinv, formMMinv ) / measure;

        for( int i = 0; i < localdim; i++ )
        for( int j = 0; j < localdim; j++ )
        {
            int index_of_entry = s * localdim * localdim + i * localdim + j;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * localdim + i;
            entry.column = s * localdim + j;
            entry.value  = fullMMinv( i, j );
            
            ret.setentry( index_of_entry, entry );
        }
        
    }
    
    return ret;
}








// inline DenseMatrix elementmassmatrix( int n, int ambientdim, int r, int k, DenseMatrix Jacobian )
// {
//     
//     assert( Jacobian.getdimin() == n && Jacobian.getdimout() == ambientdim ); 
//     assert( n <= ambientdim );
//     assert( n >= 0 && n >= k );
//         
//     DenseMatrix polyMM = polynomialmassmatrix( n, r );
//     
//     DenseMatrix formMM = SubdeterminantMatrix( mesh.getGradientProductMatrix( n, t ), k );
//         
//     return TensorProduct( polyMM, formMM ) * absolute( determinant( Jacobian ) ) / factorial_numerical( n );
//     
// }

