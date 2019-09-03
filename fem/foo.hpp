#ifndef INCLUDEGUARD_FOO
#define INCLUDEGUARD_FOO

#include <vector>
#include <algorithm>

#include "../basic.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/linearoperator.hpp"
#include "../operators/diagonaloperator.hpp"
#include "../dense/densematrix.hpp"



/************************
****
****  Class for Sparse Matrices  
****  - instantiates LinearOperator
****  
************************/


int SullivanSpanSize( int n, int k, int r )
{
    assert( 0 <= n && 0 <= k && 0 <= r );
    assert( k <= n );
    return binomial( n+1 + r, r ) * binomial( n+1, k )
}






// Generates a matrix whose columns are the barycentric coordinates
// of the interpolation points 


DenseMatrix InterpolationPointsBarycentricCoordinates( int dim, int r )
{
    const auto multi_indices = generateMultiIndices( dim, r );
    
    assert( multi_indices.size() == binomial( n+1+r, r ) );
    
    const int N = multi_indices.size();
    
    const Float delta = 0.1;
    
    DenseMatrix ret( N+1, N );
    
    for( int i = 0; i < N; i++ )
        ret.setcolumn( i, multi_indices[i].getfloatvector().shift( delta ).scaleinverse( r + (dim+1) * delta ).getentries() );

    for( int i = 0; i < N; i++ ) {
        assert( ret.getcolumn(i).isnonnegative() );
        assert( ret.getcolumn.sumnorm() > 0.999 && ret.getcolumn.sumnorm() < 1.001 );
    }
    
    return ret;
}



// Suppose that the columns of bc are the interpolation points
// in barycentric coordinates, with the coordinates in the rows.
// 
// we generate a matrix where each column represent a barycentric polynomial
// and the rows correspond to the interpolation points



DenseMatrix EvaluationMatrix( int dim, r, DenseMatrix lpsbc )
{
    const auto mis = generateMultiIndices( dim, r );
    
    assert( dim        == lpsbc.getdimout() );
    assert( mis.size() == lpsbc.getdimin()  );
    
    const int N = mis.size();
    
    DenseMatrix ret( N, N );
    
    for( int c = 0; c < N; c++ ) // c -> barycentric poly 
    for( int r = 0; r < N; r++ ) // r -> interpolation point
    {
        ret(r,c) = 1.;
        for( int d = 0; d <= dim; d++ )
            if( mis[c][d] != 0 )
                ret(r,c) *= power( lpsbc(d,r), mis[c][d] );
    }
    
    return ret;
    
}



// Given the Jacobian of the coordinate transformation, 
// this produces a matrix that transform from Euclidean coordinates 
// to barycentric coordinates 

DenseMatrix BarycentricProjectionMatrix( const DenseMatrix& Jac )
{
    assert( Jac.getdimout() >= Jac.getdimin() );
    
    // We implement J^+ = inv( J^t J ) J^t
    
    auto Jt = Tranpose(J);
    
    auto F = Inverse( Jt * J ) * Jt; // n x d
    
    DenseMatrix ret( Jac.getdimin()+1, Jac.getdimout(), 0.0 );
    
    for( int r = 0; r < Jac.getdimin();  r++ )
    for( int c = 0; c < Jac.getdimout(); c++ )
        ret( r+1, c ) = F(r,c);
    
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


FloatVector Interpolation( const Mesh& m, int dim, int k, int r, 
            std::function< FloatVector( std::vector<Float> ) > field
            )
{
    
    assert( 0 <= dim && dim <= m.getinnerdimension() );
    assert( 0 <= k && k <= dim );
    assert( 0 <= r );
    
    FloatVector ret( m.count_simplices(dim) * SullivanSpanSize() );
    
    const auto lpsbc = InterpolationPointsBarycentricCoordinates( dim, r );
    
    for( int s = 0; s < m.count_simplices(dim); s++ )
    {
        
        auto lpsbc = InterpolationPointsBarycentricCoordinates( dim, r );
        
        
        auto coords = m.getVertexCoordinateMatrix( dim, s );
        
        auto lps = coords * lpsbc;
        
        
        auto EM = EvaluationMatrix( dim, r, lpsbc );
        
        auto EMinv = Inverse( EM );
        
        
        auto Jac = m.getTransformationMatrix( dim, s );
        
        auto bpm = BarycentricProjectionMatrix( Jac );
        
        
        auto InterpolationMatrix = MatrixTensorProduct( EMinv, bpm );
        
        
        auto EvaluationVector = evaluate( lps, k, field ); // TODO: implement
        
        auto localResult = InterpolationMatrix * EvaluationVector
        
        // TODO: write local result into global vector 'ret'
        
    }
    
    return ret;
}




 

#endif
