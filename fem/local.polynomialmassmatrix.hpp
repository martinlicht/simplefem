#ifndef INCLUDEGUARD_FEM_POLYNOMIALMASSMATRIX
#define INCLUDEGUARD_FEM_POLYNOMIALMASSMATRIX


#include <vector>
#include <iostream>
// #include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/densematrix.hpp"
#include "../operators/linearoperator.hpp"




inline DenseMatrix polynomialmassmatrix( int n, int r )
{
    assert( n >= 0 && r >= 0 );
    
    // create the multiindices and the return value matrix
    
    std::vector<MultiIndex> multis = generateMultiIndices( IndexRange(0,n), r );
    
    const N = multis.size();
    
    DenseMatrix ret( N, N );
    
    // add 
    
    for( int i = 0; i < N; i++ )
    for( int j = 0; j < N; j++ )
    {
        
        MultiIndex sumindex = multis[i] + multis[j];
        
        ret( i, j ) = factorial( sumindex ) / ( n + absolute( sumindex ) );
        
    }
    
    return ret;
}







#endif
