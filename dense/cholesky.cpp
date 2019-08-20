
#include "cholesky.hpp"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <cctype>
#include <cmath>

#include "../combinatorics/generateindexmaps.hpp"
#include "densematrix.hpp"
#include "../operators/floatvector.hpp"
#include "../solver/crm.hpp"




DenseMatrix CholeskyDecomposition( const DenseMatrix& A )
{
  return CholeskyDecompositionBanachchiewicz( A );
}


DenseMatrix CholeskyDecompositionBanachchiewicz( const DenseMatrix& A )
{
  A.check();
  assert( A.issquare() );
  
  DenseMatrix L = A;
  L.set( 0. );
  const int dim = A.getdimout();
  
  for( int r = 0; r < dim; r++ ){
    
    for( int c = 0; c < r; c++ ) {
      
      L(r,c) = A(r,c);
      for( int k = 0; k < c; k++ )
        L(r,c) -= L(r,k) * L(c,k);
      L(r,c) /= L(c,c);
      
    }
    
    L(r,r) = A(r,r);
    for( int k = 0; k < r; k++ )
      L(r,r) -= L(r,k) * L(r,k);
    L(r,r) = sqrt( L(r,r) );
    
  }
  
  L.check();
  return L;
}








