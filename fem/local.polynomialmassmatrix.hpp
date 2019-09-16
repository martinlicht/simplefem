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
    
    assert( N == binomial_integer( n + r , r ) );
    
    DenseMatrix ret( N, N, notanumber );
    
    // create the entries  
    
    for( int i = 0; i < N; i++ )
    for( int j = 0; j < N; j++ )
    {
        
        MultiIndex alpha = multis[i] + multis[j];
        
        //  alpha! / (n+|alpha|)!
        assert( absolute( alpha ) == 2*r );
        
        ret( i, j ) = factorial_integer(n) * (Float)alpha.factorial() / (Float) factorial_integer( (int64_t)n + 2*r ); 

        if( ret( i, j ) <= 0. ) {
            std::cout << multis[i] << multis[j] << ret (i, j );
            unreachable();
        }
        

        //// TODO unclear whether this factor makes any sense.
	// it needs to be fixed how this matrix is scaled.
	// if it is scaled as mass matrix of unit simplex, then we later only multiple a determinant
        // if it is scaled so that the volume term is missing, then we multiply by the measure later.
        
    }
    
    assert( ret.isfinite() );
    
    //std::cout << ret << std::endl;

    if( not ret.ispositive() )
        std::cout << ret;
    
    assert( ret.isnonnegative() );
    assert( ret.ispositive() );
    
    return ret;
}







#endif
