
#include "matrixtensorproduct.hpp"

#include <cassert>



DenseMatrix MatrixTensorProduct( const DenseMatrix& A, const DenseMatrix& B )
{
    A.check(); B.check();
    int newrows = A.getdimout() * B.getdimout();
    int newcols = A.getdimin()  * B.getdimin();
    
    DenseMatrix ret( newrows, newcols );
    assert( ret.getdimout() == newrows );
    assert( ret.getdimin()  == newcols );

    const int num_A_rows = A.getdimout();
    const int num_A_cols  = A.getdimin();
    const int num_B_rows = B.getdimout();
    const int num_B_cols  = B.getdimin();
    
    for( int rl = 0; rl < num_A_rows; rl++ )
    for( int cl = 0; cl < num_A_cols; cl++ )
    for( int rr = 0; rr < num_B_rows; rr++ )
    for( int cr = 0; cr < num_B_cols; cr++ )
    {
        ret( rl * num_B_rows + rr, cl * num_B_cols + cr )
        =
        A( rl, cl ) * B( rr, cr );
    }
    
    ret.check();
    return ret;
}




