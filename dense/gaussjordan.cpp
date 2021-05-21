
#include "gaussjordan.hpp"

#include <new>

 // LR factorization, column pivot, 
 
DenseMatrix GaussJordan( DenseMatrix mat )
{
    
    assert( mat.issquare() );
    
    const int n = mat.getdimout();
//     std::vector<int> pivotcol( n, -17 );
//     std::vector<int> colperm ( n, -17 );
//     for( int c = 0; c < n; c++ ) colperm[c] = c;

    
    DenseMatrix ret(n);
    ret.unitmatrix();
    
    // 1. eliminate lower triangular part, save coeffecients
    
    for( int i = 0; i < n; i++ ) {
        
        for( int k = 0; k < n; k++ ) { // each
            
            if( i == k ) continue; 
            
            Float coeff = - mat( k, i ) / mat( i, i );
            
            for( int j = i; j < n; j++ ) {
                mat( k, j ) = mat( k, j ) + coeff * mat( i, j );
            }

            for( int j = 0; j <= i; j++ ) {
                ret( k, j ) = ret( k, j ) + coeff * ret( i, j );
            }
            
        }
        
        Float coeff = 1. / mat(i,i);
        
        for( int j = 0; j <= i; j++ ) {
            ret(i,j) *= coeff;
        }
        
        for( int j = i; j < n; j++ ) {
            mat(i,j) *= coeff;
        }
            
    }
    
// //     2. normalize the diagonals
//     
//     for( int i = 0; i < n; i++ ) {
//         
//         Float coeff = 1. / mat(i,i);
//         
//         for( int k = 0; k < n; k++ ) {
//             mat(i,k) *= coeff;
//             ret(i,k) *= coeff;
//         }
        
//         ret(i,i) = coeff;
//     }
    
    // finished!
    
    // LOG << mat;
    // LOG << ret;
    
    return ret;
}

 

 
 
 
 
DenseMatrix GaussJordanInplace( DenseMatrix mat, bool pivoting )
{
    
    assert( mat.issquare() );
    
    const int n = mat.getdimout();
    
    int* pivots = nullptr;
    if(pivoting) pivots = new (std::nothrow) int[n];
    
    for( int i = 0; i < n; i++ ) {
        
        if( pivoting ) {
            
            int c_max = i;
            for( int c = i+1; c < n; c++ )
                if( absolute(mat(i,c)) > absolute(mat(i,c_max)) ) 
                    c_max = c;
            
            pivots[i] = c_max;
            mat.swapcolumn( c_max, i );
            
        }
        
        for( int k = 0; k < n; k++ ) { // each
            
            if( i == k ) continue; 
            
            assert( absolute(mat(i,i)) != 0.0 );
            
            Float coeff = - mat( k, i ) / mat( i, i );
            
            for( int j = i+1; j < n; j++ ) {
                mat( k, j ) = mat( k, j ) + coeff * mat( i, j );
            }

            mat( k, i ) = coeff;
            
            for( int j = 0; j < i; j++ ) {
                mat( k, j ) = mat( k, j ) + coeff * mat( i, j );
            }
            
        }
        
        Float coeff = 1. / mat(i,i);
        
        for( int j = 0; j < i; j++ ) {
            mat(i,j) *= coeff;
        }
        
        mat(i,i) = coeff;
        
        for( int j = i+1; j < n; j++ ) {
            mat(i,j) *= coeff;
        }
            
    }
    
    if( pivoting ) {
        for( int i = n-1; i >= 0; i-- )
//         for( int i = 0; i < n; i++ ) 
        {
//             LOG << "swap " << i << space << pivots[i] << nl;
            mat.swaprow( i, pivots[i] );
        }
    }
    
    // finished!
    
    if( pivoting ) delete[] pivots;
    
    return mat;
}











