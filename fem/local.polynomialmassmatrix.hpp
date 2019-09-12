#ifndef INCLUDEGUARD_FEM_POLYNOMIALMASSMATRIX
#define INCLUDEGUARD_FEM_POLYNOMIALMASSMATRIX


#include <vector>
#include <iostream>
// #include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"
#include "../operators/linearoperator.hpp"




inline DenseMatrix polynomialmassmatrix( int n, int r )
{
    assert( n >= 0 && r >= 0 );
    
    // create the multiindices and the return value matrix
    
    std::vector<MultiIndex> multis = generateMultiIndices( IndexRange(0,n), r );
    
    const int N = multis.size();
    
    assert( N == binomial( n + r , r ) );
    
    DenseMatrix ret( N, N );
    
    // create the entries  
    
    for( int i = 0; i < N; i++ )
    for( int j = 0; j < N; j++ )
    {
        
        MultiIndex alpha = multis[i] + multis[j];
        
        //  alpha! / (n+|alpha|)!
        
        ret( i, j ) = factorial(n) * alpha.factorial() / (Float) factorial( n + absolute( alpha ) );
        
    }
    
    
    //std::cout << ret << std::endl;
        
    
    return ret;
}







#endif
