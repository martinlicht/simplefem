
#include "matrixtensorproduct.hpp"

#include <cctype>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <list>




DenseMatrix MatrixTensorProduct( const DenseMatrix& left, const DenseMatrix& right )
{
    left.check(); right.check();
    int newrows = left.getdimout() * right.getdimout();
    int newcols = left.getdimin() * right.getdimin();
    
    DenseMatrix ret( newrows, newcols );
    assert( ret.getdimout() == newrows );
    assert( ret.getdimin() == newcols );
    
    for( int rl = 0; rl < left.getdimout(); rl++ )
    for( int cl = 0; cl < left.getdimin(); cl++ )
    for( int rr = 0; rr < right.getdimout(); rr++ )
    for( int cr = 0; cr < right.getdimin(); cr++ )
    {
        ret( rl * right.getdimout() + rr, cl * right.getdimin() + cr )
        =
        left( rl, cl ) * right( rr, cr );
    }
    
    ret.check();
    return ret;
}




