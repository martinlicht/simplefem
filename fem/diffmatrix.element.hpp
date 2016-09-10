#ifndef INCLUDEGUARD_DIFFMATRIX_ELEMENT
#define INCLUDEGUARD_DIFFMATRIX_ELEMENT


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
****  Method produces an element diff matrix 
****  
****  - workds completely stand alone 
****  - any dimension 
****  - 
****  
****  
*******************/



DenseMatrix calculateElementDiffMatrix(
                int, // innerdimension
                int, // polynomial degree
                int // form degree 
                );

inline int polysigmaindex2fullindex( int n, int r, int k, int poly_i, int sigma_i )
{
  int poly_n = binomial( n + r, r );
  int sigma_n = binomial( n + 1, k );
  assert( 0 <= poly_i && poly_i < poly_n );
  assert( 0 <= sigma_i && sigma_i < sigma_n );
  return poly_i * poly_n + sigma_i;
}

inline int fullindex2polysigmaindex( int n, int r, int k, int full_i, int& poly_i, int& sigma_i )
{
  int poly_n = binomial( n + r, r );
  int sigma_n = binomial( n + 1, k );
  assert( 0 <= full_i && full_i < poly_n * sigma_n );
  poly_i = full_i / poly_n;
  sigma_i = full_i % poly_n;
}





#endif