#ifndef INCLUDEGUARD_MASSMATRIX_ELEMENT
#define INCLUDEGUARD_MASSMATRIX_ELEMENT


#include <vector>
#include <iostream>
// #include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/densematrix.hpp"
#include "../operators/linearoperator.hpp"




/*******************
****  
****  Method produces an element mass matrix 
****  
****  - workds completely stand alone 
****  - any dimension 
****  - 
****  
****  
*******************/


// Calculates the mass matrix for the barycentric polynomials 
Float integrateBarycentricPolynomialUnitsimplex( MultiIndex alpha );



// Calculates the mass matrix for the barycentric polynomials 
DenseMatrix calculateScalarMassMatrixUnitSimplex(
                int // polynomial degree 
                );

// Calculates the mass matrix for the barycentric exterior derivatives 
std::vector<FloatVector> calculateBarycentricDiffs(
                int, int, // inner and outer dimension 
                const std::vector<FloatVector>& // vertices of simplex 
                );

// Calculates the mass matrix for the barycentric exterior derivatives 
DenseMatrix calculateBDproductMatrix(
                int, int, // inner and outer dimension 
                const std::vector<FloatVector>&, // vertices of simplex 
                int // form degree 
                );

DenseMatrix calculateElementMassMatrix(
                int, int, // inner and outer dimension 
                const std::vector<FloatVector>&, // vertices of simplex 
                int, // polynomial degree 
                int // form degree 
                );

Float simplexvolume( int, int, const std::vector<FloatVector>& );


std::vector<FloatVector> calculateSimplexheightpoints(
                int innerdim, int outerdim, 
                const std::vector<FloatVector>& vertices 
                );
                
std::vector<FloatVector> calculateSimplexheightvectors(
                int innerdim, int outerdim, 
                const std::vector<FloatVector>& vertices 
                );
                




#endif