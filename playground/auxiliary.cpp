

// #include <cassert>
#include <iostream>
#include <vector>
#include <iterator>

#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/matrixalgorithm.hpp"
#include "auxiliary.hpp"




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




std::vector<FloatVector> calculateBarycentricDiffs(
                int innerdim, int outerdim, 
                const std::vector<FloatVector>& vertices 
                )
{
    assert( outerdim >= 0 && innerdim >= 0 ); 
    assert( vertices.size() + 1 == innerdim );
    
    std::vector<FloatVector> ret = calculateSimplexheightvectors(innerdim, outerdim, vertices);
    
    for( int v = 0; v < innerdim+1; v++ ) 
        ret[v] /= std::pow( ret[v].norm(), 2. );
    
    return ret;
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
    
    DenseMatrix BD( outerdim, innerdim+1, calculateBarycentricDiffs( innerdim, outerdim, vertices ) );
    DenseMatrix BDSP = BD.transpose() * BD;
    return SubdeterminantMatrix( BDSP, formdegree );
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
    DenseMatrix form_part = calculateBDproductMatrix( innerdim, outerdim, vertices, formdegree );
    Float volume = simplexvolume( innerdim, outerdim, vertices );
    
    return volume * MatrixTensorProduct( poly_part, form_part );
    
}



Float simplexvolume( int innerdim, int outerdim, const std::vector<FloatVector>& vertices )
{
    
    assert( outerdim >= 0 && innerdim >= 0 ); 
    assert( vertices.size() + 1 == innerdim );
    
    /* Volume formula according to 
    http://www.mathpages.com/home/kmath664/kmath664.htm
    */
    
    DenseMatrix temp( outerdim, innerdim+1, vertices );
    temp = temp.transpose() * temp;
    temp.add( 1. );
    Float fakinnerdim = factorial( innerdim );
    return temp.determinant() / ( fakinnerdim * fakinnerdim );
}




std::vector<FloatVector> calculateSimplexheightpoints(
                int innerdim, int outerdim, 
                const std::vector<FloatVector>& vertices 
                )
{
    assert( outerdim >= 0 && innerdim >= 0 ); 
    assert( vertices.size() + 1 == innerdim );
    
    std::vector<FloatVector> ret( innerdim+1, FloatVector(outerdim) );

    /* TODO */
    /* https://www.math.auckland.ac.nz/~waldron/Preprints/Barycentric/barycentric.pdf */
    
    return ret;
}


std::vector<FloatVector> calculateSimplexheightvectors(
                int innerdim, int outerdim, 
                const std::vector<FloatVector>& vertices 
                )
{
    assert( outerdim >= 0 && innerdim >= 0 ); 
    assert( vertices.size() + 1 == innerdim );
    
    std::vector<FloatVector> ret = calculateSimplexheightpoints( innerdim, outerdim, vertices );
    
    for( int v = 0; v < innerdim+1; v++ ) 
        ret[v] = vertices[v] - ret[v];
    
    return ret;
}

