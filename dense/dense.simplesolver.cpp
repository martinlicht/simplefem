
#include "dense.simplesolver.hpp"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <cctype>

#include "../combinatorics/generateindexmaps.hpp"
#include "densematrix.hpp"
#include "../operators/floatvector.hpp"
#include "../solver/crm.hpp"





DenseMatrix DiagonalPart( const DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    DenseMatrix Ret( A.getdimout(), A.getdimin(), 0. );
    for( int c = 0; c < A.getdimin(); c++ )
        Ret(c,c) = A(c,c); 
    return Ret;
}

DenseMatrix DiagonalInverse( const DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    DenseMatrix Ret( A.getdimout(), A.getdimin(), 0. );
    for( int c = 0; c < A.getdimin(); c++ )
        Ret(c,c) = 1. / A(c,c); 
    return Ret;
}

void InvertDiagonal( DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    for( int c = 0; c < A.getdimin(); c++ )
        A(c,c) = 1. / A(c,c); 
}

Float DiagonalDeterminant( const DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    Float ret = 1.;
    for( int c = 0; c < A.getdimin(); c++ )
        ret = ret * A(c,c);
    return ret;
}

void DiagonalSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b )
{
    A.check();
    assert( A.issquare() );
    assert( x.getdimension() == A.getdimin() );
    assert( b.getdimension() == A.getdimout() );
    for( int c = 0; c < A.getdimin(); c++ )
        x[c] = b[c] / A(c,c); 
}









DenseMatrix LowerTriangularPart( const DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    DenseMatrix Ret( A.getdimout(), A.getdimin(), 0. );
    for( int r = 0; r < A.getdimout(); r++ )
        for( int c = 0; c <= r; c++ )
            Ret(r,c) = A(r,c); 
    return Ret;
}

DenseMatrix LowerUnitTriangularPart( const DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    DenseMatrix Ret( A.getdimout(), A.getdimin(), 0. );
    for( int r = 0; r < A.getdimout(); r++ ) {
        Ret(r,r) = 1.;
        for( int c = 0; c < r; c++ )
            Ret(r,c) = A(r,c); 
    }
    return Ret;
}

DenseMatrix UpperTriangularPart( const DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    DenseMatrix Ret( A.getdimout(), A.getdimin(), 0. );
    for( int r = 0; r < A.getdimout(); r++ )
        for( int c = r; c < A.getdimin(); c++ )
            Ret(r,c) = A(r,c); 
    return Ret;
}

DenseMatrix UpperUnitTriangularPart( const DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    DenseMatrix Ret( A.getdimout(), A.getdimin(), 0. );
    for( int r = 0; r < A.getdimout(); r++ ) {
        Ret(r,r) = 1.;
        for( int c = r+1; c < A.getdimin(); c++ )
            Ret(r,c) = A(r,c); 
    }
    return Ret;
}








DenseMatrix LowerTriangularInverse( const DenseMatrix& A )
{
    assert( A.issquare() );
    DenseMatrix Ret( A );
    Ret.zeromatrix();
    
    for( int r = 0; r < A.getdimout(); r++ ) {
        for( int c = 0; c < r; c++ ) {
            Ret(r,c) = 0.;
            for( int k = c; k < r; k++ )
                Ret(r,c) = Ret(r,c) - A(r,k) * Ret(k,c);
            Ret(r,c) /= A(r,r);
        }
        Ret(r,r) = 1. / A(r,r);
    }
    Ret.check();
    return Ret;
}

DenseMatrix LowerUnitTriangularInverse( const DenseMatrix& A, bool writediagonalones )
{
    assert( A.issquare() );
    DenseMatrix Ret( A );
    Ret.zeromatrix();
    
    for( int r = 0; r < A.getdimout(); r++ ) {
        for( int c = 0; c < r; c++ ) {
            Ret(r,c) = 0.;
            for( int k = c; k < r; k++ )
                Ret(r,c) = Ret(r,c) - A(r,k) * Ret(k,c);
        }
        if( writediagonalones ) Ret(r,r) = 1.;
    }
    Ret.check();
    return Ret;
}

DenseMatrix UpperTriangularInverse( const DenseMatrix& A )
{
    assert( A.issquare() );
    DenseMatrix Ret( A );
    Ret.zeromatrix();
    
    for( int r = A.getdimout()-1; r >= 0; r-- ) {
        for( int c = A.getdimout()-1; c > r; c-- ) {
            Ret(r,c) = 0.;
            for( int k = r+1; k < c; k++ )
                Ret(r,c) -= A(r,k) * Ret(k,c);
            Ret(r,c) /= A(r,r);
        }
        Ret(r,r) = 1. / A(r,r);
    }
    
    Ret.check();
    return Ret;
}

DenseMatrix UpperUnitTriangularInverse( const DenseMatrix& A, bool writediagonalones )
{
    assert( A.issquare() );
    DenseMatrix Ret( A );
    Ret.zeromatrix();
    
    for( int r = A.getdimout()-1; r >= 0; r-- ) {
        for( int c = A.getdimout()-1; c > r; c-- ) {
            Ret(r,c) = 0.;
            for( int k = r+1; k < c; k++ )
                Ret(r,c) -= A(r,k) * Ret(k,c);
        }
        if( writediagonalones ) Ret(r,r) = 1.;
    }
    
    Ret.check();
    return Ret;
}













void InvertLowerTriangular( DenseMatrix& A )
{
    assert( A.issquare() );
    for( int r = 0; r < A.getdimout(); r++ ) {
        for( int c = 0; c < r; c++ ) {
            Float temp = 0.;
            for( int k = c; k < r; k++ )
                    temp -= A(r,k) * A(k,c);
            A(r,c) = temp / A(r,r);
        }
        A(r,r) = 1. / A(r,r);
    }
}

void InvertLowerUnitTriangular( DenseMatrix& A, bool writediagonalones )
{
    assert( A.issquare() );
    for( int r = 0; r < A.getdimout(); r++ ) {
        for( int c = 0; c < r; c++ ) {
            Float temp = 0.;
            for( int k = c; k < r; k++ )
                    temp -= A(r,k) * A(k,c);
            A(r,c) = temp;
        }
        if( writediagonalones ) A(r,r) = 1.;
    }
}

void InvertUpperTriangular( DenseMatrix& A )
{
    assert( A.issquare() );
    for( int r = A.getdimout()-1; r >= 0; r-- ) {
        for( int c = A.getdimout()-1; c > r; c-- ) {
            Float temp = 0.;
            for( int k = r+1; k < c; k++ )
                    temp -= A(r,k) * A(k,c);
            A(r,c) = temp / A(r,r);
        }
        A(r,r) = 1. / A(r,r);
    }
}

void InvertUpperUnitTriangular( DenseMatrix& A, bool writediagonalones )
{
    assert( A.issquare() );
    for( int r = A.getdimout()-1; r >= 0; r-- ) {
        for( int c = A.getdimout()-1; c > r; c-- ) {
            Float temp = 0.;
            for( int k = r+1; k < c; k++ )
                    temp -= A(r,k) * A(k,c);
            A(r,c) = temp;
        }
        if( writediagonalones ) A(r,r) = 1.;
    }
}










Float LowerTriangularDeterminant( const DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    Float ret = 1.;
    for( int c = 0; c < A.getdimin(); c++ )
        ret = ret * A(c,c);
    return ret;
}

Float UpperTriangularDeterminant( const DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    Float ret = 1.;
    for( int c = 0; c < A.getdimin(); c++ )
        ret = ret * A(c,c);
    return ret;
}










void LowerTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b )
{
    A.check();
    x.check();
    b.check();
    assert( A.issquare() );
    assert( x.getdimension() == A.getdimin() );
    assert( b.getdimension() == A.getdimout() );
    for( int r = 0; r < A.getdimout(); r++ ) {
        x[r] = b[r];
        for( int c = 0; c < r-1; c++ )
            x[r] = x[r] - A(r,c) * x[c];
        x[r] = x[r] / A(r,r);
    }
}

void LowerUnitTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b )
{
    A.check();
    x.check();
    b.check();
    assert( A.issquare() );
    assert( x.getdimension() == A.getdimin() );
    assert( b.getdimension() == A.getdimout() );
    for( int r = 0; r < A.getdimout(); r++ ) {
        x[r] = b[r];
        for( int c = 0; c < r-1; c++ )
            x[r] = x[r] - A(r,c) * x[c];
    }
}

void UpperTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b )
{
    A.check();
    x.check();
    b.check();
    assert( A.issquare() );
    assert( x.getdimension() == A.getdimin() );
    assert( b.getdimension() == A.getdimout() );
    for( int r = A.getdimout()-1; r >= 0; r-- ) {
        x[r] = b[r];
        for( int c = r+1; c < A.getdimin(); c++ )
            x[r] = x[r] - A(r,c) * x[c];
        x[r] = x[r] / A(r,r);
    }
}

void UpperUnitTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b )
{
    A.check();
    x.check();
    b.check();
    assert( A.issquare() );
    assert( x.getdimension() == A.getdimin() );
    assert( b.getdimension() == A.getdimout() );
    for( int r = A.getdimout()-1; r >= 0; r-- ) {
        x[r] = b[r];
        for( int c = r+1; c < A.getdimin(); c++ )
            x[r] = x[r] - A(r,c) * x[c];
    }
}











