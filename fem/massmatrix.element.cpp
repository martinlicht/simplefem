

#include <cassert>
#include <iostream>
#include <vector>
#include <iterator>

#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "massmatrix.element.hpp"




Float integrateBarycentricPolynomialUnitsimplex( MultiIndex alpha )
{
    int dimension = alpha.getIndexRange().cardinality();
    return factorial(dimension) * alpha.factorial() 
            / factorial( dimension + alpha.absolute() );
}




DenseMatrix calculateScalarMassMatrixUnitSimplex( int dimension, int polydegree )
{
    assert( polydegree >= 0 );
    IndexRange ir( 0, dimension );
    std::vector<MultiIndex> mis = generateMultiIndices( ir, polydegree );
    DenseMatrix ret( mis.size() );
    for( int l = 0; l < mis.size(); l++ ) 
        for( int r = 0; r < mis.size(); r++ )
            ret( l, r ) = integrateBarycentricPolynomialUnitsimplex( mis[l] + mis[r] );
    return ret;        
}




DenseMatrix calculateBarycentricDiffs(
                int innerdim, int outerdim, 
                const std::vector<FloatVector>& vertices 
                )
{
    assert( outerdim >= 0 && innerdim >= 0 ); 
    assert( vertices.size() + 1 == innerdim );
    
    DenseMatrix ret( outerdim, innerdim+1 );
    /* The differentials d\lambda_0, \dots, d\lambda_m of the 
       barycentric coordinates \lambda_0, \dots, \lambda_m. */
    /* TODO: Either with Polar decomposition or solving a linear system */
    /* Probably Polar decomposition is the better choice */
}
                

                
                
DenseMatrix calculateBDproductMatrix(
                int innerdim, int outerdim,
                const std::vector<FloatVector>& vertices,
                int formdegree
                )
{
    assert( outerdim >= 0 && innerdim >= 0 ); 
    assert( vertices.size() + 1 == innerdim );
    assert( 0 <= formdegree && formdegree <= innerdim );
    
    DenseMatrix BD = calculateBarycentricDiffs( innerdim, outerdim, vertices );
    DenseMatrix BDSP = BD.transpose() * BD;
    return Subdeterminantmatrix( BDSP, formdegree );
}




DenseMatrix calculateElementMassMatrix(
                int innerdim, int outerdim,
                const std::vector<FloatVector>& vertices,
                int polydegree,
                int formdegree 
                )
{
    assert( outerdim >= 0 && innerdim >= 0 ); 
    assert( vertices.size() + 1 == innerdim );
    assert( 0 <= formdegree && formdegree <= innerdim );
    assert( 0 <= polydegree && polydegree <= innerdim );
    
    DenseMatrix poly_part = calculateScalarMassMatrixUnitSimplex( innerdim, polydegree );
    DenseMatrix form_part = calculateScalarMassMatrixUnitSimplex( innerdim, outerdim, vertices, formdegree );
    Float volume = simplexvolume( vertices );
    
    return volume * MatrixTensorProduct( poly_part, form_part );
    
}
