
#include "dense.functions.hpp"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <cctype>

#include "../combinatorics/generateindexmaps.hpp"
#include "densematrix.hpp"
#include "qr.factorization.hpp"
#include "lu.factorization.hpp"
#include "../operators/floatvector.hpp"
#include "../solver/crm.hpp"





DenseMatrix Transpose( const DenseMatrix& src ) 
{
    src.check();
    DenseMatrix ret( src.getdimin(), src.getdimout() );
    for( int r = 0; r < src.getdimout(); r++ )
    for( int c = 0; c < src.getdimin(); c++ )
        ret(c,r) = ret(r,c);
    ret.check();
    return ret;
}

DenseMatrix TransposeSquare( const DenseMatrix& src ) 
{
    src.check();
    return Transpose( src );
}

void TransposeInSitu( DenseMatrix& src )
{
  src.check();
  const int numrows = src.getdimout();
  const int numcols = src.getdimin();
  
  for( unsigned int start = 0; start < numcols * numrows; ++start )
  {
    
    unsigned int next = start;
    unsigned int i = 0;
    
    do {
      ++i;
      next = (next % numrows) * numcols + next / numrows;
    } while (next > start);

    if ( next >= start && i != 1 )
    {
      const double tmp = src( start / numcols, start % numcols );
      next = start;
      do {
        i = (next % numrows) * numcols + next / numrows;
        src( next / numcols, next % numcols ) = ( i == start ) ? tmp : src( i / numcols, i % numcols );
        next = i;
      } while (next > start);
    }
  
  }
  
}

void TransposeSquareInSitu( DenseMatrix& src ) 
{
    src.check();
    for( int r = 0; r < src.getdimout(); r++ )
    for( int c = r+1; c < src.getdimin(); c++ )
        { Float t = src(c,r); src(c,r) = src(r,c); src(r,c) = t; };
    src.check();
}

//         ret.entries.at( c * src.getdimout() + r ) = src.entries.at( r * src.getdimin() + c );




Float MatrixTrace( const DenseMatrix& src )
{
    src.check();
    assert( src.issquare() );
    Float ret = 0.;
    for( int i = 0; i < src.getdimout(); i++ )
      ret += src(i,i);
    return ret;
}




DenseMatrix Gerschgorin( const DenseMatrix& src )
{
    src.check();
    assert( src.issquare() );
    return GerschgorinRow( src );
}

DenseMatrix GerschgorinRow( const DenseMatrix& src )
{
    src.check();
    assert( src.issquare() );
    DenseMatrix ret( src.getdimout(), 2 );
    for( int r = 0; r < src.getdimout(); r++ )
    {
        ret( r, 0 ) = src(r,r);
        ret( r, 1 ) = 0.;
        for( int c = 0; c < r; c++ )
            ret( r, 1 ) += absolute( ret(r,c) );
        for( int c = r+1; c < src.getdimout(); c++ )
            ret( r, 1 ) += absolute( ret(r,c) );
    }
    return ret;
}

DenseMatrix GerschgorinColumn( const DenseMatrix& src )
{
    src.check();
    assert( src.issquare() );
    DenseMatrix ret( src.getdimout(), 2 );
    for( int c = 0; c < src.getdimout(); c++ )
    {
        ret( c, 0 ) = src(c,c);
        ret( c, 1 ) = 0.;
        for( int r = 0; r < c; r++ )
            ret( c, 1 ) += absolute( ret(r,c) );
        for( int r = c+1; r < src.getdimout(); r++ )
            ret( c, 1 ) += absolute( ret(r,c) );
    }
    return ret;
}






Float NormL1( const DenseMatrix& src )
{
    src.check();
    Float ret = 0.;
    for( int r = 0; r < src.getdimout(); r++ )
    for( int c = 0; c < src.getdimin(); c++ )
        ret += absolute( src(r,c) );
    return ret;
}

Float NormFrobenius( const DenseMatrix& src )
{
    src.check();
    Float ret = 0.;
    for( int r = 0; r < src.getdimout(); r++ )
    for( int c = 0; c < src.getdimin(); c++ )
        ret += absolute( src(r,c) ) * absolute( src(r,c) );
    ret = sqrt( ret );
    return ret;
}

Float NormMax( const DenseMatrix& src )
{
    src.check();
    Float ret = 0.;
    for( int r = 0; r < src.getdimout(); r++ )
    for( int c = 0; c < src.getdimin(); c++ )
        ret = maximum( ret, absolute( src(r,c) ) );
    return ret;
}

Float NormLp( const DenseMatrix& src, Float p )
{
    src.check();
    assert( 1. <= p );
    Float ret = 0.;
    for( int r = 0; r < src.getdimout(); r++ )
    for( int c = 0; c < src.getdimin(); c++ )
        ret += pow( absolute( src(r,c) ), p );
    ret = pow( ret, 1. / p );
    return ret;
}

Float NormRowCol( const DenseMatrix& src, Float p, Float q )
{
    src.check();
    assert( 1. <= p && 1. <= q );
    Float ret = 0.;
    for( int r = 0; r < src.getdimout(); r++ ) {
        Float zeile = 0.;
        for( int c = 0; c < src.getdimin(); c++ )
            zeile += pow( absolute( src(r,c) ), q );
        ret += pow( zeile, p/q );
    }
    ret = pow( ret, 1. / p );
    return ret;
}

Float NormColRow( const DenseMatrix& src, Float p, Float q )
{
    src.check();
    assert( 1. <= p && 1. <= q );
    Float ret = 0.;
    for( int c = 0; c < src.getdimin(); c++ ) {
        Float spalte = 0.;
        for( int r = 0; r < src.getdimout(); r++ )
            spalte += pow( absolute( src(r,c) ), q );
        ret += pow( spalte, p/q );
    }
    ret = pow( ret, 1. / p );
    return ret;
}

Float NormOperatorL1( const DenseMatrix& src )
{
    src.check();
    Float ret = 0.;
    for( int c = 0; c < src.getdimin(); c++ ) {
        Float spalte = 0.;
        for( int r = 0; r < src.getdimout(); r++ )
            spalte += absolute( src(r,c) );
        ret = maximum( ret, spalte );
    }
    return ret;
}

Float NormOperatorMax( const DenseMatrix& src )
{
    src.check();
    Float ret = 0.;
    for( int r = 0; r < src.getdimout(); r++ ) {
        Float zeile = 0.;
        for( int c = 0; c < src.getdimin(); c++ )
            zeile += absolute( src(r,c) );
        ret = maximum( ret, zeile );
    }
    return ret;
}








void InverseAndDeterminant( const DenseMatrix& A, DenseMatrix& Ainv, Float& Adet )
{
    DenseMatrix Q( A ), R( A );
    QRFactorization( A, Q, R );
    Ainv = UpperTriangularInverse(R) * Transpose(Q);
    Adet = UpperTriangularDeterminant(R);
}


DenseMatrix Inverse( const DenseMatrix& A )
{
    DenseMatrix Q( A ), R( A );
    QRFactorization( A, Q, R );
    return UpperTriangularInverse(R) * Transpose(Q);
}

Float Determinant( const DenseMatrix& A )
{
    DenseMatrix Q( A );
    DenseMatrix R( A );
    QRFactorization( A, Q, R );
    return UpperTriangularDeterminant( R );
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




