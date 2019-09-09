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
        
        for( int k = 0; k < n; k++ ) {
            
            if( i == k ) continue; 
            
            Float coeff = - mat( k, i ) / mat( i, i );
            
            for( int j = 0; j < n; j++ ) {
                mat( k, j ) = mat( k, j ) + coeff * mat( i, j );
                ret( k, j ) = ret( k, j ) + coeff * ret( i, j );
            }
            
//             mat( k, i ) = coeff;
            
        }
            
    }
    
    // 2. normalize the diagonals
    
    for( int i = 0; i < n; i++ ) {
        
        Float coeff = 1. / mat(i,i);
        
        for( int k = 0; k < n; k++ ) {
            mat(i,k) *= coeff;
            ret(i,k) *= coeff;
        }
        
    }
    
    // finished!
    
    std::cout << mat;
    std::cout << ret;
    
    return ret;
}

 


#endif
