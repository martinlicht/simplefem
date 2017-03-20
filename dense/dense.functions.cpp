
#include "dense.functions.hpp"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <cctype>

#include "../combinatorics/generateindexmaps.hpp"
#include "densematrix.hpp"
#include "dense.factorization.hpp"
#include "../operators/floatvector.hpp"
#include "../solver/crm.hpp"



void InverseAndDeterminant( const DenseMatrix& A, DenseMatrix& Ainv, Float& Adet )
{
    DenseMatrix Q( A ), R( A );
    QRFactorization( A, Q, R );
    Ainv = UpperTriangularInverse(R) * Q.transpose();
    Adet = UpperTriangularDeterminant(R);
}


DenseMatrix Inverse( const DenseMatrix& A )
{
    DenseMatrix Q( A ), R( A );
    QRFactorization( A, Q, R );
    return UpperTriangularInverse(R) * Q.transpose();
}

Float Determinant( const DenseMatrix& A )
{
    DenseMatrix Q( A );
    DenseMatrix R( A );
    QRFactorization( A, Q, R );
    return UpperTriangularDeterminant( R );
}



Float TriangularDeterminant( const DenseMatrix& src )
{
    assert( src.issquare() );
    Float ret = 1.;
    for( int i = 0; i < src.getdimout(); i++ )
        ret *= src(i,i);
    return ret;
}



DenseMatrix UpperTriangularInverse( const DenseMatrix& src )
{
    assert( src.issquare() );
    DenseMatrix ret( src );
    ret.zeromatrix();
    const int D = src.getdimout();
    for( int c = D-1; c >= 0; c-- ) {
        ret(c,c) = 1. / src(c,c);
        for( int r = c-1; r >= 0; r-- ) {
            ret(r,c) = 0.;
            for( int k = c; k > r; k-- )
                    ret(r,c) -= src(r,k) * ret(k,c);
            ret(r,c) /= src(r,r);
        }
    }
    ret.check();
    return ret;
}


Float UpperTriangularDeterminant( const DenseMatrix& A )
{
    assert( A.issquare() );
    Float ret = 1.;
    for( int c = 0; c < A.getdimin(); c++ )
        ret *= A(c,c);
    return ret;
}


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



DenseMatrix Subdeterminantmatrix( const DenseMatrix& A, int k )
{
    A.check();
    assert( A.issquare() );
    assert( 0 <= k && k <= A.getdimin() );
    
    const int n = A.getdimin();
    IndexRange fromrange = IndexRange( 0, k-1 );
    IndexRange torange = IndexRange( 0, n-1 );
    std::vector<IndexMap> sigmas = generateSigmas( fromrange, torange );
    
    DenseMatrix ret( sigmas.size() );
    for( int rim = 0; rim < sigmas.size(); rim++ )
    for( int cim = 0; cim < sigmas.size(); cim++ )
    {
        ret(rim,cim) = determinant( A.submatrix( sigmas.at(rim), sigmas.at(cim) ) );
    }
    
    ret.check();
    return ret;
}




