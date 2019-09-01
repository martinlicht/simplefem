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
#include "../dense/densematrix.hpp"
#include "../dense/matrixtensorproduct.hpp"
#include "../operators/linearoperator.hpp"

#include "../fem/local.polynomialmassmatrix.hpp"







inline SparseMatrix FEECBrokenMassMatrix( Mesh& mesh, int n, int k, int r )
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
    
    for( int s = 0; s < num_simplices; s++ )
    {
        
        DenseMatrix Jac = mesh.getTransformationJacobian( n, s );
            
        DenseMatrix formMM = SubdeterminantMatrix( mesh.getGradientProductMatrix( n, s ), k );
    
        DenseMatrix fullMM = MatrixTensorProduct( polyMM, formMM ) * absolute( determinant( Jac ) ) / factorial( n );
        
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
