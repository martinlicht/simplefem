#ifndef INCLUDEGUARD_FOO
#define INCLUDEGUARD_FOO

#include <vector>
#include <algorithm>

#include "../basic.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/linearoperator.hpp"
#include "../operators/diagonaloperator.hpp"
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






// Generates a matrix whose columns are the barycentric coordinates
// of the interpolation points 

inline DenseMatrix InterpolationPointsBarycentricCoordinates( int n, int r )
{
    assert( 0 <= n && 0 <= r );
    
    const auto multi_indices = generateMultiIndices( IndexRange(0,n), r );
    
    assert( multi_indices.size() == binomial_integer( n+r, r ) );
    
    const Float delta = +0.10000000;
    
    DenseMatrix ret( n+1, multi_indices.size(), 0. );
    
    for( int i = 0; i < multi_indices.size(); i++ )
        ret.setcolumn( i, FloatVector( multi_indices[i].getvalues() ).shift( delta ).scaleinverse( r + (n+1) * delta ) );

    if( r != 0 )
    for( int i = 0; i < multi_indices.size(); i++ ) {
        assert( ret.getcolumn(i).isnonnegative() );
        assert( ret.getcolumn(i).sumnorm() > 0.9999 && ret.getcolumn(i).sumnorm() < 1.0001 );
    }
    
    return ret;
}



// Suppose that the columns of lpsbc are the interpolation points
// in barycentric coordinates.
// 
// we generate a matrix where each column represent a barycentric polynomial
// and the rows correspond to the interpolation points

inline DenseMatrix EvaluationMatrix( int dim, int r, const DenseMatrix& lpsbc )
{
    assert( 0 <= dim and 0 <= r );
    
    const auto mis = generateMultiIndices( IndexRange(0,dim), r );
    
    assert( dim+1      == lpsbc.getdimout() );
    assert( mis.size() == lpsbc.getdimin()  );
    
    const int N = mis.size();
    
    DenseMatrix ret( N, N );
    
    for( int c = 0; c < N; c++ ) // c -> barycentric poly 
    for( int r = 0; r < N; r++ ) // r -> interpolation point
    {
        ret(r,c) = 1.;
        for( int d = 0; d <= dim; d++ )
            if( mis[c][d] != 0 )
                ret(r,c) *= power( lpsbc(d,r), (Float) mis[c][d] );
    }
    
    return ret;
    
}



// Suppose that the columns of lpsbc are the interpolation points
// in barycentric coordinates.
// 
// we generate a matrix where each column represent a barycentric polynomial
// and the rows correspond to the interpolation points

// TODO Reread and compare with previous version 

inline DenseMatrix EvaluationMatrix( std::vector<MultiIndex> mis, const DenseMatrix& lpsbc )
{
    
    const int num_points = lpsbc.getdimin();
    const int num_mis    = mis.size();
    
    int dim = lpsbc.getdimout() - 1;
        
    DenseMatrix ret( num_points, num_mis );
    
    for( int c = 0; c < num_mis;    c++ ) // c -> barycentric poly 
    for( int r = 0; r < num_points; r++ ) // r -> interpolation point
    {
        assert( mis[c].getSourceRange().max() == dim+1 );
        assert( mis[c].getSourceRange().min() == 0     );

        ret(r,c) = 1.;
        for( int d = 0; d <= dim; d++ )
            if( mis[c][d] != 0 )
                ret(r,c) *= power( lpsbc(d,r), (Float) mis[c][d] );
    }
    
    return ret;
    
}



// Given the Jacobian of the coordinate transformation, 
// this produces a matrix that transform from Euclidean coordinates 
// to barycentric coordinates 

inline DenseMatrix BarycentricProjectionMatrix( const DenseMatrix& J )
{
    assert( J.getdimout() >= J.getdimin() );
    
    // We implement J^+ = inv( J^t J ) J^t
    
    auto Jt = Transpose(J);
    
    auto F = Inverse( Jt * J ) * Jt; // n x d
    
    DenseMatrix ret( J.getdimin()+1, J.getdimout(), 0.0 );
    
    for( int r = 0; r < J.getdimin();  r++ )
    for( int c = 0; c < J.getdimout(); c++ )
        ret( r+1, c ) = F(r,c);
    
    return ret;
}


inline DenseMatrix BarycentricProjectionMatrixNEUERVERSUCH( const DenseMatrix& J )
{
    assert( J.getdimout() >= J.getdimin() );
    
    auto F = Transpose(J);
    
    DenseMatrix ret( J.getdimin()+1, J.getdimout(), 0.0 );
    
    for( int r = 0; r < J.getdimin();  r++ )
    for( int c = 0; c < J.getdimout(); c++ )
        ret( r+1, c ) = F(r,c);
    
    return ret;
}





inline DenseMatrix EvaluateField( 
            int dim, int k, int r, 
            const DenseMatrix& lps, 
            std::function< FloatVector( const FloatVector& ) > field
            )
{
    
    assert( 0 <= dim );
    assert( 0 <= k && k <= dim );
    assert( 0 <= r );
    assert( lps.getdimout() == dim );
    
    const auto fielddim = binomial_integer(dim,k);
    
    DenseMatrix ret( fielddim, lps.getdimin() );
    
    for( int p = 0; p < lps.getdimin(); p++ )
    {
        auto point = lps.getcolumn(p);
        
        auto value = field( point );
        
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
    
    //const auto degree_r_mis = generateMultiIndices( IndexRange(0,dim), r );
    
    const auto lpsbc = InterpolationPointsBarycentricCoordinates( dim, r );
    
    const auto EM = EvaluationMatrix( dim, r, lpsbc );
        
    const auto EMinv = Inverse( EM );
    
    assert( EM.isfinite()    );
    assert( EMinv.isfinite() );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < m.count_simplices(dim); s++ )
    {
        
        const auto coords = m.getVertexCoordinateMatrix( dim, s );
        
        const auto lps = coords * lpsbc;
        
//         std::cout << lps << space << lpsbc << coords << nl;
        
        assert( lps.isfinite() );
        
        const auto Jac = m.getTransformationJacobian( dim, s );
        
        const auto bpm = BarycentricProjectionMatrixNEUERVERSUCH( Jac );
        
        const auto P = SubdeterminantMatrix( bpm, k ); //MatrixAlternatingPower( bpm, k );
        
        assert( P.isfinite()    );
        
        const auto InterpolationMatrix = MatrixTensorProduct( EMinv, P );
        assert( InterpolationMatrix.isfinite()    );
        
        const auto Evaluations = EvaluateField( dim, k, r, lps, field );
        const auto EvaluationVector = Evaluations.flattencolumns();
        
        assert( Evaluations.isfinite()      );
        assert( EvaluationVector.isfinite() );
        

//         std::cout << InterpolationMatrix.getdimin() << space << EvaluationVector.getdimension() << nl;
        
        const auto localResult = InterpolationMatrix * EvaluationVector;

        assert( localResult.isfinite()    );
        // std::cout << EvaluationVector << std::endl;

        if( k == 0 ) {
        
            assert( ( EM * localResult - EvaluationVector ).issmall() );
            assert( ( InterpolationMatrix - EMinv ).issmall() );
        }
        
        assert( localResult.getdimension() == SullivanSpanSize(dim,k,r) );
        
        ret.setslice( s * SullivanSpanSize(dim,k,r), localResult );
        
    }
    
    return ret;
}




 

#endif
