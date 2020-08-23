#ifndef INCLUDEGUARD_FOO
#define INCLUDEGUARD_FOO

#include <algorithm>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/linearoperator.hpp"
#include "../operators/simpleoperators.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/matrixtensorproduct.hpp"
#include "../dense/functions.hpp"
#include "../mesh/mesh.hpp"


/************************
****
****  Class for Sparse Matrices  
****  - instantiates LinearOperator
****  
************************/


inline int SullivanSpanSize( int n, int k, int r )  __attribute__ ((const));

inline int SullivanSpanSize( int n, int k, int r )
{
    assert( 0 <= n && 0 <= k && 0 <= r );
    assert( k <= n );
    return binomial_integer( n + r, r ) * binomial_integer( n+1, k );
}






// 
// Generates a matrix whose columns are the barycentric coordinates
// of the interpolation points over an n dimensional simplex 
// such that polynomials of degree r can be reconstructed.
// 
// The barycentric coordinates are chosen such that 
// the points are located in the interior of the simplex.
// (See parameter delta)
// 
// Size of returned matrix:
// [n+1] x [ n+r choose r ]

inline DenseMatrix InterpolationPointsBarycentricCoordinates( int n, int r )
{
    assert( 0 <= n && 0 <= r );
    
    const auto multi_indices = generateMultiIndices( IndexRange(0,n), r );
    
    assert( multi_indices.size() == binomial_integer( n+r, r ) );
    
    const Float delta = +0.10000000;
    
    DenseMatrix ret( n+1, multi_indices.size(), 0. );
    
    assert( delta > 0 );
    
    for( int i = 0; i < multi_indices.size(); i++ )
        ret.setcolumn( i, FloatVector( multi_indices[i].getvalues() ).shift( delta ).scaleinverse( r + (n+1) * delta ) );

    if( r != 0 )
    for( int i = 0; i < multi_indices.size(); i++ ) {
        assert( ret.getcolumn(i).isnonnegative() );
        assert( ret.getcolumn(i).sumnorm() > 0.9999 && ret.getcolumn(i).sumnorm() < 1.0001 );
    }
    
    return ret;
}



// 
// Suppose that the columns of lpsbc are the interpolation points
// (in barycentric coordinates) to interpolate degree r polynomials 
// over a d simplex exactly.
// 
// The output matrix is as follows:
// - rows correspond to the interpolation points 
// - columns correspond to the multiindices (Lagrange basis)
// - the entries are value of the corresponding polynomial at the corresponding point
// 
// Size of returned matrix:
// [ n+r choose r ] x [ n+r choose r ]

inline DenseMatrix EvaluationMatrix( int dim, int r, const DenseMatrix& lpsbc )
{
    assert( 0 <= dim and 0 <= r );
    
    const auto mis = generateMultiIndices( IndexRange(0,dim), r );
    
    assert( dim+1      == lpsbc.getdimout()            );
    assert( mis.size() == lpsbc.getdimin()             );
    assert( mis.size() == binomial_integer( dim+r, r ) );
    const int N = mis.size();
    
    DenseMatrix ret( N, N );
    
    for( int col = 0; col < N; col++ ) // col -> barycentric poly 
    for( int row = 0; row < N; row++ ) // row -> interpolation point 
    {
        ret(row,col) = 1.;
        for( int d = 0; d <= dim; d++ )
            if( mis[col][d] != 0 )
                ret(row,col) *= power_numerical( lpsbc(d,row), mis[col][d] );
    }
    
    return ret;
    
}



// 
// Suppose that the columns of lpsbc are the interpolation points
// (in barycentric coordinates) to interpolate degree r polynomials 
// over a d simplex exactly.
// 
// Suppose that mis cointains a collection of multiindices
// over a d simplex 
// 
// The output matrix is as follows:
// - columns correspond to the interpolation points 
// - rows correspond to the multiindices (Lagrange basis)
// - the entries are value of the corresponding polynomial at the corresponding point
// 
// TODO Reread and compare with previous version 

inline DenseMatrix EvaluationMatrix( std::vector<MultiIndex> mis, const DenseMatrix& lpsbc )
{
    
    const int num_points = lpsbc.getdimin();
    const int num_mis    = mis.size();    
    const int dim        = lpsbc.getdimout() - 1;
        
    DenseMatrix ret( num_points, num_mis );
    
    for( int c = 0; c < num_mis;    c++ ) // c -> barycentric poly 
    for( int r = 0; r < num_points; r++ ) // r -> interpolation point
    {
        assert( mis[c].getSourceRange().max() == dim+1 );
        assert( mis[c].getSourceRange().min() == 0     );

        ret(r,c) = 1.;
        for( int d = 0; d <= dim; d++ )
            if( mis[c][d] != 0 )
                ret(r,c) *= power_numerical( lpsbc(d,r), mis[c][d] );
    }
    
    return ret;
    
}



// Given the Jacobian of the coordinate transformation, 
// this produces a matrix that transform from Euclidean coordinates 
// to barycentric coordinates 

// inline DenseMatrix BarycentricProjectionMatrixALTERVERSUCH( const DenseMatrix& J )
// {
//     assert( J.getdimout() >= J.getdimin() );
    
//     // We implement J^+ = inv( J^t J ) J^t
    
//     const auto Jt = Transpose(J);
    
//     const auto F = Inverse( Jt * J ) * Jt; // n x d
    
//     DenseMatrix ret( J.getdimin()+1, J.getdimout(), 0.0 );
    
//     for( int r = 0; r < J.getdimin();  r++ )
//     for( int c = 0; c < J.getdimout(); c++ )
//         ret( r+1, c ) = F(r,c);
    
//     return ret;
// }













// 
// The input matrix J is a Jacobian of the transformation which transforms 
// the [dimin] reference simplex onto a simplex 
// within [dimout] dimensional ambient space 
// 
// We transpose the matrix and prepend a zero row.
// 
// If we multiply a vector in standard basis coordinates 
// we get the same vector in barycentric coordinates, 
// using only the gradients 1 through n. 
// Gradient 0 is redundant and thus the zero-th row is 0.
// 
// Size of output matrix:
// [ J.dimin()+1 x J.dimout() ]
// 


inline DenseMatrix BarycentricProjectionMatrix( const DenseMatrix& J )
{
    assert( J.getdimout() >= J.getdimin() );
    
    const auto F = Transpose(J);
    
    DenseMatrix ret( J.getdimin()+1, J.getdimout(), 0.0 );
    
    for( int r = 0; r < J.getdimin();  r++ )
    for( int c = 0; c < J.getdimout(); c++ )
        ret( r+1, c ) = F(r,c);
    
    return ret;
}






// 
// evaluate a field at given physical points 
// and collect its values 
// 
// Given a k-differential formin dim dimensional ambient space 
// and lps being a matrix of size [outerdimension] x [var1]
// we let the output be as follows.
// 
// The output is a matrix of size [dim choose k] x [var1]
// whose column are the evaluations of the field 
// at each of the interpolation points.
// 

inline DenseMatrix EvaluateField( 
            int dim, int k,  
            const DenseMatrix& lps, 
            std::function< FloatVector( const FloatVector& ) > field
            )
{
    
    assert( 0 <= dim );
    assert( 0 <= k && k <= dim );
    assert( lps.getdimout() == dim );
    
    const auto fielddim = binomial_integer(dim,k);
    
    DenseMatrix ret( fielddim, lps.getdimin() );
    
    for( int p = 0; p < lps.getdimin(); p++ )
    {
        const auto point = lps.getcolumn(p);
        
        const auto value = field( point );
        
        assert( value.getdimension() == binomial_integer( dim, k ) );
        
        ret.setcolumn( p, value );
    }
    
    return ret;
    
}













// TODO:
// 
// Given a function in ambient space, giving forms in Euclidean coordinates
// 
// 1. Evaluate the form at the Lagrange points 
// 
// 2. Get the coefficients of the barycentric polynomials
// 
// 3. Transform the Euclidean components to their barycentric projections
// 
// The output is a vector whose coefficients represent the volume-wise interpolations 
// into the local Sullivan space of k-forms with polynomial degree r 
// 


inline FloatVector Interpolation( 
            const Mesh& m, 
            int dim, int k, int r, 
            std::function< FloatVector( const FloatVector& ) > field
            )
{
    
    assert( 0 <= dim && dim <= m.getinnerdimension() );
    assert( 0 <= k && k <= dim );
    assert( 0 <= r );
    
    FloatVector ret( m.count_simplices(dim) * SullivanSpanSize(dim,k,r) );
    
    
    // lpsbc are the barycentric coordinates to interpolate 
    // degree r polynomials over a dim simplex 
    // [dim+1] x [number of lagrange points] 
    // 
    // EM is the values of the Lagrange polynomials at those points 
    // [n+r choose r] x [number of lagrange points] (actually, square)
    // 
    // EMinv is the inverse of EM.
    // [same size as EM]
    
    const auto lpsbc = InterpolationPointsBarycentricCoordinates( dim, r );
    
    const auto EM = EvaluationMatrix( dim, r, lpsbc );
        
    const auto EMinv = Inverse( EM );
    
    assert( lpsbc.getdimin()  == SullivanSpanSize(dim,0,r) );
    assert( lpsbc.getdimout() == dim+1 );
    assert( EM.getdimout()    == SullivanSpanSize(dim,0,r) );
    assert( EM.getdimin()     == SullivanSpanSize(dim,0,r) );
    assert( EM.isfinite()    );
    assert( EMinv.isfinite() );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < m.count_simplices(dim); s++ )
    {
        
        // We obtain the coordinate matrix (outerdim) x (dim+1) 
        // of the dim-dimensional simplex with index s.
        // Each column contains the physical coordinates of the vertices
        // 
        // lps are the physical coordinates of the interpolation points 
        // (outerdim) x (number of Lagrange points)
        // 
        
        const auto coords = m.getVertexCoordinateMatrix( dim, s );
        
        const auto lps = coords * lpsbc;
        
        // 
        // Jac is a matrix of size [outerdimension] x [dim]
        // which is the jacobian of the transformation mapping of simplex s
        // 
        // bpm is a matrix of size [dim+1]x[outerdimension]
        // that is the barycentric projection matrix (see paper)
        // 
        // P is the matrix of subdeterminants of size k 
        // of that matrix 
        // 
        
        const auto Jac = m.getTransformationJacobian( dim, s );
        
        const auto bpm = BarycentricProjectionMatrix( Jac );
        
        const auto P = SubdeterminantMatrix( bpm, k );
        
        // The InterpolationMatrix is the tensor product matrix 
        // of EMinv and P.
        // 
        // Evaluations is a matrix of size [dim choose k] x [lps.dimin]
        // whose columsn are the outputs of the field at the Lagrange points
        // 
        // localResults contains the coeffecients of the interpolation
        // in the canonical basis (Lagrange polynomials x barycentric gradients)
        
        const auto InterpolationMatrix = MatrixTensorProduct( EMinv, P );
        
        const auto Evaluations      = EvaluateField( dim, k, lps, field );
        const auto EvaluationVector = Evaluations.flattencolumns();
        
        const auto localResult = InterpolationMatrix * EvaluationVector;

        #ifndef NDEBUG

        assert( lps.isfinite() );
        
        assert( P.isfinite()    );
        
        assert( InterpolationMatrix.isfinite()    );

        assert( Evaluations.isfinite()      );
        assert( EvaluationVector.isfinite() );
        
        assert( localResult.isfinite()    );

        if( k == 0 ) {
            if( !( EM * localResult - EvaluationVector ).issmall() ) {
                LOG << EM * localResult << space << EvaluationVector << nl;
                LOG << EM * localResult - EvaluationVector << nl;
                LOG << ( EM * localResult - EvaluationVector ).norm() << nl;
            }
            assert( ( EM * localResult - EvaluationVector ).issmall() );
            assert( ( InterpolationMatrix - EMinv ).issmall() );
        }
        #endif
        
        assert( localResult.getdimension() == SullivanSpanSize(dim,k,r) );
        
        ret.setslice( s * SullivanSpanSize(dim,k,r), localResult );
        
    }
    
    return ret;
}




 

#endif
