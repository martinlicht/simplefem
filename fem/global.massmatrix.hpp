#ifndef INCLUDEGUARD_FEM_FEECBROKENMASSMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENMASSMATRIX


#include <vector>
#include <iostream>
// #include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/cholesky.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/matrixtensorproduct.hpp"
#include "../operators/linearoperator.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/local.polynomialmassmatrix.hpp"







inline SparseMatrix FEECBrokenMassMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim = binomial( n+r, n ) * binomial( n, k );
    
    const int dim_in  = num_simplices * localdim;
    const int dim_out = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    DenseMatrix polyMM = polynomialmassmatrix( n, r );
    
//     std::cout << polyMM << std::endl;
        
    for( int s = 0; s < num_simplices; s++ )
    {
        
        Float measure = mesh.getMeasure( n, s );
            
        DenseMatrix formMM = SubdeterminantMatrix( mesh.getGradientProductMatrix( n, s ), k );
    
        DenseMatrix fullMM = MatrixTensorProduct( polyMM, formMM ) * measure;
        
//         std::cout << measure << std::endl;
        
//         std::cout << formMM << std::endl;
        
//         std::cout << fullMM << std::endl;
        
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





inline SparseMatrix FEECBrokenMassMatrixRightFactor( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim = binomial( n+r, n ) * binomial( n, k );
    
    const int dim_in  = num_simplices * localdim;
    const int dim_out = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    DenseMatrix polyMM = polynomialmassmatrix( n, r );
    
    DenseMatrix polyMM_right = Transpose(CholeskyDecomposition(polyMM));
    
//     std::cout << polyMM << std::endl;
        
    for( int s = 0; s < num_simplices; s++ )
    {
        
        Float measure = mesh.getMeasure( n, s );
            
        DenseMatrix formMM_right = SubdeterminantMatrix( mesh.getGradientProductMatrixRightFactor( n, s ), k );
    
        DenseMatrix fullMM_right = MatrixTensorProduct( polyMM_right, formMM_right ) * sqrt(measure);
        
        {
            
            /* TEST WHETHER THE ARITHMETICS WORK OUT */
            
            DenseMatrix formMM = SubdeterminantMatrix( mesh.getGradientProductMatrix( n, s ), k );
            DenseMatrix fullMM = MatrixTensorProduct( polyMM, formMM ) * measure;
            
            assert( ( Transpose(polyMM_right) * polyMM_right - polyMM ).issmall() ); 
            assert( ( Transpose(formMM_right) * formMM_right - formMM ).issmall() ); 
            assert( ( Transpose(fullMM_right) * fullMM_right - fullMM ).issmall() ); 
            
        }
        
        
        for( int i = 0; i < localdim; i++ )
        for( int j = 0; j < localdim; j++ )
        {
            int index_of_entry = s * localdim * localdim + i * localdim + j;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * localdim + i;
            entry.column = s * localdim + j;
            entry.value  = fullMM_right( i, j );
            
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
//     return TensorProduct( polyMM, formMM ) * absolute( determinant( Jacobian ) ) / factorial( n );
//     
// }



#endif
