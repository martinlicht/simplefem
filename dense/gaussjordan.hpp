#ifndef INCLUDEGUARD_DENSE_LU_FACTORIZATION
#define INCLUDEGUARD_DENSE_LU_FACTORIZATION

#include <vector>

#include "../basic.hpp"

#include "densematrix.hpp"


 
 // LR factorization, column pivot, 
 
inline DenseMatrix GaussJordan( DenseMatrix mat )
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
    
    std::cout << mat;
    std::cout << ret;
    
    return ret;
}

 

 
 
 
 
inline DenseMatrix GaussJordanInplace( DenseMatrix mat )
{
    
    assert( mat.issquare() );
    
    const int n = mat.getdimout();
    
    for( int i = 0; i < n; i++ ) {
        
        for( int k = 0; k < n; k++ ) { // each
            
            if( i == k ) continue; 
            
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
    
    
    // finished!
    
    return mat;
}

 

#endif
