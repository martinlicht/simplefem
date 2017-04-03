
#include "cholesky.hpp"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <cctype>

#include "../combinatorics/generateindexmaps.hpp"
#include "densematrix.hpp"
#include "../operators/floatvector.hpp"
#include "../solver/crm.hpp"





DenseMatrix CholeskyDecomposition( const DenseMatrix& src )
{
    src.check();
    assert( src.issquare() );
    
    DenseMatrix ret = src;
    ret.set( 0. );
    const int D = src.getdimout();
    
    for( int k = 0; k < D; k++ ){
        ret(k,k) = src(k,k);
        for( int j = 0; j < k; j++ )
            ret(k,k) -= ret(j,k) * ret(j,k);
        ret(k,k) = sqrt( ret(k,k) );
        for( int i = k+1; i < D; i++ ) {
            ret(k,i) = src(k,i);
            for( int j = 0; j < k; j++ )
                    ret(k,i) -= ret(j,i) * ret(j,k);
            ret(k,i) /= ret(k,k);
        }
    }
    
    ret.check();
    return ret;
}






