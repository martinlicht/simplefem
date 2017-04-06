
#include "iterativeinverse.hpp"
#include "scalarfunctions.hpp"

void newtoniteration( DenseMatrix& X, const DenseMatrix& A, int n )
{
  assert( n >= 0 );
  assert( A.issquare() );
  assert( X.issquare() );
  assert( X.getdimout() == A.getdimout() );
  
  Float alpha = EigenvalueEstimate( A );
  
  X = ( 1. / ( alpha * alpha ) ) * A;
  
  for( int i = 0; i < n; i++ )
    X = ( 2 * X ) - ( X * ( A * X ) );
}
