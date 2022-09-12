#ifndef INCLUDEGUARD_FEM_UTILITIES_HPP
#define INCLUDEGUARD_FEM_UTILITIES_HPP


#include "../basic.hpp"
#include "../dense/densematrix.hpp"
#include "../mesh/mesh.hpp"





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

DenseMatrix InterpolationPointsBarycentricCoordinates( int n, int r );



// 
// Suppose that the columns of bcs are evaluation points
// (in barycentric coordinates) over a d simplex.
// 
// The output matrix is as follows:
// - row indices correspond to the evaluation points 
// - column indices correspond to the multiindices (standard Lagrange basis of degree r)
// - the entries are value of the corresponding polynomial at the corresponding point
// 
// Size of returned matrix:
// [ n+r choose r ] x [ number of evaluation points ]
// 
// NOTE: if EM is the output matrix and x some coefficient vector for a polynomial,
//       then y = EM * x is the values of that polynomial at the Lagrange points
//       If EM is invertible, then x = inv(EM) y provides the monomial coefficients 
//       for a given vector y of point values.

DenseMatrix EvaluationMatrix( int r, const DenseMatrix& bcs );







// 
// Suppose that the columns of bcs are evaluation points
// (in barycentric coordinates) over a d simplex.
// 
// The output matrix is as follows:
// - row indices correspond to the multiindices (standard Lagrange basis of degree r)
// - column indices correspond to the evaluation points 
// - each entry is the coefficient of the corresponding multindex 
//   for the Lagrange polynomial associated to the evaluation point 
//   TODO: fix whether that should be transposed.
// 
// Size of returned matrix:
// [ n+r choose r ] x [ n+r choose r ]

DenseMatrix LagrangePolynomialCoefficients( int n, int r );




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
// TODO: This should be a member of the mesh class rather than a function 
//       in the FEM module. Should be moved from here.

DenseMatrix BarycentricProjectionMatrix( const DenseMatrix& J );




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
// at each of the evaluation points, given in Euclidean coordinates.
// 

DenseMatrix EvaluateField( 
            int outerdim, int k,  
            const DenseMatrix& lps, 
            std::function< FloatVector( const FloatVector& ) > field
            );





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

FloatVector Interpolation( 
            const Mesh& m, 
            int dim, int k, int r, 
            std::function< FloatVector( const FloatVector& ) > field
            );


#endif
