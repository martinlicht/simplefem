
#include "dense.scalarfunctions.hpp"

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
            ret( r, 1 ) += absolute( src(r,c) );
        for( int c = r+1; c < src.getdimout(); c++ )
            ret( r, 1 ) += absolute( src(r,c) );
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
            ret( c, 1 ) += absolute( src(r,c) );
        for( int r = c+1; r < src.getdimout(); r++ )
            ret( c, 1 ) += absolute( src(r,c) );
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
    for( int c = 0; c < src.getdimin(); c++ ) {
        Float av = absolute( src(r,c) );
        ret += av * av;
    }
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














