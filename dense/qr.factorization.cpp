
#include "qr.factorization.hpp"

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

void LQFactorization( const DenseMatrix& A, DenseMatrix& L, DenseMatrix& Q )
{
    A.check();
    L.check();
    Q.check();
    
    // TODO: check algorithm
    assert(false);
    
    assert( A.getdimout() == L.getdimout() );
    assert( A.getdimin()  == Q.getdimin()  ); 
    assert( Q.getdimout() == L.getdimin() );
    assert( A.getdimout() <= A.getdimin() );
    assert( L.issquare() );
    
    L.zeromatrix();
    
    for( int r = 0; r < A.getdimout(); r++ ) {
        FloatVector v = A.getrow(r);
        for( int i = 0; i < r; i++ ){
                L(r,i) = v * Q.getrow(i);
                v -= L(r,i) * Q.getrow(i);
        }
        L(r,r) = sqrt( v*v );
        Q.setrow( r, v / L(r,r) );
    }
    
}


// // // void QRFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t )
// // // {
// // //     if( t == 0 )
// // //         return;
// // //     if( t == 1 )
// // //         QRFactorization( A, Q, R );
// // //     else {
// // //         DenseMatrix Qw(Q), Qv(Q);
// // //         DenseMatrix Rw(R), Rv(R);
// // //         QRFactorizationRepeated( A, Qw, Rw, t-1 );
// // //         QRFactorization( Qw, Qv, Rv );
// // //         Q = Qv;
// // //         R = Rv * Rw;
// // //     }
// // // }
// // // 
// // // 
// // // void LQFactorizationRepeated( const DenseMatrix& A, DenseMatrix& L, DenseMatrix& Q, unsigned int t )
// // // {
// // //     assert(false);
// // //     
// // //     if( t == 0 )
// // //         return;
// // //     if( t == 1 )
// // //         LQFactorization( A, L, Q );
// // //     else {
// // //         DenseMatrix Qw(Q), Qv(Q);
// // //         DenseMatrix Lw(L), Lv(L);
// // //         LQFactorizationRepeated( A, Lw, Qw, t-1 );
// // //         LQFactorization( Qw, Lv, Qv );
// // //         Q = Qv;
// // //         L = Lw * Lv;
// // //     }
// // // }





