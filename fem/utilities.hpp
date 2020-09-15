#ifndef INCLUDEGUARD_FEM_UTILITIES_HPP
#define INCLUDEGUARD_FEM_UTILITIES_HPP

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
    for( const auto& mi : multi_indices ) assert( mi.absolute() == r );
    
    const Float delta = 0.1;
    
    DenseMatrix ret( n+1, static_cast<int>( multi_indices.size() ) );
    
    assert( ret.getdimout() == n+1 );
    assert( ret.getdimin() == multi_indices.size() );
    
    assert( delta > 0 );
    
    for( int i = 0; i < multi_indices.size(); i++ )
        ret.setcolumn( i, FloatVector( multi_indices[i].getvalues() ).shift( delta ).scaleinverse( r + (n+1) * delta ) );

    for( int i = 0; i < ret.getdimout(); i++ ) {
        assert( ret.getcolumn(i).isnonnegative() );
        assert( ret.getcolumn(i).sumnorm() > 0.9999 && ret.getcolumn(i).sumnorm() < 1.0001 );
    }
    
    assert( ret.isfinite() );
    
    return ret;
}

















// 
// Suppose that the columns of bcs are evaluation points
// (in barycentric coordinates) over a d simplex.
// 
// The output matrix is as follows:
// - rows correspond to the evaluation points 
// - columns correspond to the multiindices (standard Lagrange basis of degree r)
// - the entries are value of the corresponding polynomial at the corresponding point
// 
// Size of returned matrix:
// [ n+r choose r ] x [ number of evaluation points ]

inline DenseMatrix EvaluationMatrix( int r, const DenseMatrix& bcs )
{
    assert( 0 <= r );
    assert( 0 < bcs.getdimout() );
    
    const int dim = bcs.getdimout()-1;
    
    const auto mis = generateMultiIndices( IndexRange(0,dim), r );
    
    const int number_of_evaluation_points = bcs.getdimin();
    
    const int number_of_polynomials = mis.size();
    
    assert( mis.size() == binomial_integer( dim+r, r ) );
    
    DenseMatrix ret( number_of_polynomials, number_of_evaluation_points );
    
    for( int row = 0; row < number_of_evaluation_points; row++ ) // row -> interpolation point 
    for( int col = 0; col < number_of_polynomials;       col++ ) // col -> barycentric poly 
    {
        ret(row,col) = 1.;
        for( int d = 0; d <= dim; d++ )
            if( mis[col][d] != 0 )
                ret(row,col) *= power_numerical( bcs(d,row), mis[col][d] );
    }
    
    assert( ret.isfinite() );
    
    return ret;
    
}





























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
    
    assert( ret.isfinite() );
    
    return ret;
}























// 
// evaluate a field at given physical points 
// and collect its values 
// 
// Given a k-differential form in dim dimensional ambient space 
// and lps being a matrix of size [outerdimension] x [var1]
// we let the output be as follows.
// 
// Here, [var1] is the number of evaluation points.
// 
// The output is a matrix of size [dim choose k] x [var1]
// whose column are the evaluations of the field 
// at each of the evaluation points.
// 

inline DenseMatrix EvaluateField( 
            int outerdim, int k,  
            const DenseMatrix& lps, 
            std::function< FloatVector( const FloatVector& ) > field
            )
{
    
    assert( 0 <= outerdim );
    assert( 0 <= k && k <= outerdim );
    assert( lps.getdimout() == outerdim );
    
    const auto fielddim = binomial_integer(outerdim,k);
    
    const auto number_of_evaluation_points = lps.getdimin();
    
    DenseMatrix ret( fielddim, number_of_evaluation_points );
    
    for( int p = 0; p < number_of_evaluation_points; p++ )
    {
        const auto evaluation_point = lps.getcolumn(p);
        
        const auto value = field( evaluation_point );
        
        assert( value.getdimension() == binomial_integer( outerdim, k ) );
        
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
    
    const int outerdim = m.getouterdimension();
    
    // lpsbc are the barycentric coordinates to interpolate 
    // degree r polynomials over a dim simplex 
    // [dim+1] x [number of lagrange points] 
    // 
    // EM is the function values of the Lagrange polynomials at those points 
    // [n+r choose r] x [number of lagrange points] (actually, square)
    // 
    // EMinv is the inverse of EM.
    // [same size as EM]
    
    const auto lpsbc = InterpolationPointsBarycentricCoordinates( dim, r );
    
    const auto EM = EvaluationMatrix( r, lpsbc );
        
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
        
        const auto Evaluations      = EvaluateField( outerdim, k, lps, field );
        const auto EvaluationVector = Evaluations.flattencolumns();
        
        const auto localResult = InterpolationMatrix * EvaluationVector;

        #ifndef NDEBUG

        assert( lps.isfinite() );
        
        assert( P.isfinite() );
        
        assert( InterpolationMatrix.isfinite() );

        assert( Evaluations.isfinite()      );
        assert( EvaluationVector.isfinite() );
        
        assert( localResult.isfinite() );

        if( k == 0 ) {
            
            assert( ( P - DenseMatrix( 1, 1, 1. ) ).iszero() );
            assert( ( InterpolationMatrix - EMinv ).iszero() );
            
            if( not ( EM * localResult - EvaluationVector ).issmall() ) {
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
