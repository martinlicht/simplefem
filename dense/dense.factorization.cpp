
#include "dense.factorization.hpp"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <cctype>

#include "../combinatorics/generateindexmaps.hpp"
#include "densematrix.hpp"
#include "../operators/floatvector.hpp"
#include "../solver/crm.hpp"



void QRFactorization( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R )
{
    A.check();
    Q.check();
    R.check();
    assert( A.getdimout() == Q.getdimout() );
    assert( Q.getdimin() == R.getdimout() );
    assert( A.getdimin() == R.getdimin() ); 
    assert( A.getdimin() <= A.getdimout() );
    assert( R.issquare() );
    
    R.zeromatrix();
    
    for( int c = 0; c < A.getdimin(); c++ ) {
        FloatVector u = A.getcolumn(c);
        for( int j = 0; j < c; j++ ){
                R(j,c) = u * Q.getcolumn(j);
                u -= R(j,c) * Q.getcolumn(j);
        }
        R(c,c) = sqrt( u*u );
        Q.setcolumn( c, u / R(c,c) );
    }
    
}


void QRFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t )
{
    if( t == 0 )
        return;
    if( t == 1 )
        QRFactorization( A, Q, R );
    else {
        DenseMatrix Qw(Q), Qv(Q);
        DenseMatrix Rw(R), Rv(R);
        QRFactorizationRepeated( A, Qw, Rw, t-1 );
        QRFactorization( Qw, Qv, Rv );
        Q = Qv;
        R = Rv * Rw;
    }
}




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






