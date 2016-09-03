#ifndef INCLUDEGUARD_MASSMATRIX_ELEMENT
#define INCLUDEGUARD_MASSMATRIX_ELEMENT


#include <vector>
#include <iostream>
#include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
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
DenseMatrix calculateScalarMassMatrix(
                int, int, // inner and outer dimension 
                vector<floatvector>, // vertices of simplex 
                int, // polynomial degree 
                );

// Calculates the mass matrix for the barycentric exterior derivatives 
DenseMatrix calculateVectorMassMatrix(
                int, int, // inner and outer dimension 
                vector<floatvector> // vertices of simplex 
                );

DenseMatrix calculateElementMassMatrix(
                int, int, // inner and outer dimension 
                vector<floatvector>, // vertices of simplex 
                int, // polynomial degree 
                int // form degree 
                );



#endif