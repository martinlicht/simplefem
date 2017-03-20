
#include "dense.simplesolver.hpp"

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

DenseMatrix InvertDiagonal( DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    DenseMatrix Ret( A.getdimout(), A.getdimin(), 0. );
    for( int c = 0; c < A.getdimin(); c++ )
        A(c,c) = 1. / A(c,c); 
    return Ret;
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

// DenseMatrix LowerTriangularInverse( const DenseMatrix& A )
// {
//     assert( false ); // TODO: Change order here. 
//     assert( A.issquare() );
//     DenseMatrix ret( A );
//     ret.zeromatrix();
//     const int D = A.getdimout();
//     for( int c = D-1; c >= 0; c-- ) {
//         ret(c,c) = 1. / A(c,c);
//         for( int r = c-1; r >= 0; r-- ) {
//             ret(r,c) = 0.;
//             for( int k = c; k > r; k-- )
//                     ret(r,c) -= A(r,k) * ret(k,c);
//             ret(r,c) /= A(r,r);
//         }
//     }
//     ret.check();
//     return ret;
// }

Float LowerTriangularDeterminant( const DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );
    Float ret = 1.;
    for( int c = 0; c < A.getdimin(); c++ )
        ret = ret * A(c,c);
    return ret;
}

// void LowerTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );





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

DenseMatrix UpperTriangularInverse( const DenseMatrix& A )
{
    assert( A.issquare() );
    DenseMatrix Ret( A );
    Ret.zeromatrix();
    
    const int D = A.getdimout();
    
    for( int c = D-1; c >= 0; c-- ) {
        Ret(c,c) = 1. / A(c,c);
        for( int r = c-1; r >= 0; r-- ) {
            Ret(r,c) = 0.;
            for( int k = c; k > r; k-- )
                    Ret(r,c) -= A(r,k) * Ret(k,c);
            Ret(r,c) /= A(r,r);
        }
    }
    
    Ret.check();
    return Ret;
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

// void UpperTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );








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

// DenseMatrix LowerUnitTriangularInverse( const DenseMatrix& A );

// void LowerUnitTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );






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

// DenseMatrix UpperUnitTriangularInverse( const DenseMatrix& A );

// void UpperUnitTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );



// // // // // DenseMatrix UpperTriangularInverse( const DenseMatrix& A )
// // // // // {
// // // // //     assert( A.issquare() );
// // // // //     DenseMatrix ret( A );
// // // // //     ret.zeromatrix();
// // // // //     const int D = A.getdimout();
// // // // //     for( int c = D-1; c >= 0; c-- ) {
// // // // //         ret(c,c) = 1. / A(c,c);
// // // // //         for( int r = c-1; r >= 0; r-- ) {
// // // // //             ret(r,c) = 0.;
// // // // //             for( int k = c; k > r; k-- )
// // // // //                     ret(r,c) -= A(r,k) * ret(k,c);
// // // // //             ret(r,c) /= A(r,r);
// // // // //         }
// // // // //     }
// // // // //     ret.check();
// // // // //     return ret;
// // // // // }





