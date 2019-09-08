
#include "functions.hpp"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <cctype>

#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/heappermgen.hpp"
#include "densematrix.hpp"
#include "qr.factorization.hpp"
#include "lu.factorization.hpp"
#include "../operators/floatvector.hpp"
#include "../solver/crm.hpp"





DenseMatrix Transpose( const DenseMatrix& src ) 
{
    src.check();
    DenseMatrix ret( src.getdimin(), src.getdimout() );
    for( int r = 0; r < src.getdimout(); r++ )
    for( int c = 0; c < src.getdimin(); c++ )
        ret(c,r) = src(r,c);
    ret.check();
    return ret;
}

DenseMatrix TransposeSquare( const DenseMatrix& src ) 
{
    src.check();
    return Transpose( src );
}

void TransposeInSitu( DenseMatrix& src )
{
  src.check();
  const int numrows = src.getdimout();
  const int numcols = src.getdimin();
  
  for( unsigned int start = 0; start < numcols * numrows; ++start )
  {
    
    unsigned int next = start;
    unsigned int i = 0;
    
    do {
      ++i;
      next = (next % numrows) * numcols + next / numrows;
    } while (next > start);

    if ( next >= start && i != 1 )
    {
      const double tmp = src( start / numcols, start % numcols );
      next = start;
      do {
        i = (next % numrows) * numcols + next / numrows;
        src( next / numcols, next % numcols ) = ( i == start ) ? tmp : src( i / numcols, i % numcols );
        next = i;
      } while (next > start);
    }
  
  }
  
}

void TransposeSquareInSitu( DenseMatrix& src ) 
{
    src.check();
    for( int r = 0; r < src.getdimout(); r++ )
    for( int c = r+1; c < src.getdimin(); c++ )
        { Float t = src(c,r); src(c,r) = src(r,c); src(r,c) = t; };
    src.check();
}

//         ret.entries.at( c * src.getdimout() + r ) = src.entries.at( r * src.getdimin() + c );











Float Determinant( const DenseMatrix& A )
{
    assert( A.issquare() );
    
    if( A.getdimin() == 0 ) {
        
        return 1.;
        
    } else if( false && A.getdimin() == 1 ) {
        
        return A(0,0);
        
    } else if( false && A.getdimin() == 2 ) {
        
        return A(0,0) * A(1,1) - A(0,1) * A(1,0);
        
    } else if( false && A.getdimin() == 3 ) {
        
        return + A(0,0) * A(1,1) * A(2,2) // 1 2 3 + 
               - A(0,0) * A(1,2) * A(2,1) // 1 3 2 - 
               - A(0,1) * A(1,0) * A(2,2) // 2 1 3 - 
               + A(0,1) * A(1,2) * A(2,0) // 2 3 1 + 
               - A(0,2) * A(1,1) * A(2,0) // 3 2 1 - 
               + A(0,2) * A(1,0) * A(2,1) // 3 1 2 + 
               ;
        
    } else if( false && A.getdimin() == 4 ) {
        
        return + A(0,0) * A(1,1) * A(2,2) * A(3,3) // 0 1 2 3 + 
               - A(0,0) * A(1,2) * A(2,1) * A(3,3) // 0 2 1 3 - 
               - A(0,1) * A(1,0) * A(2,2) * A(3,3) // 1 0 2 3 - 
               + A(0,1) * A(1,2) * A(2,0) * A(3,3) // 1 2 0 3 + 
               - A(0,2) * A(1,1) * A(2,0) * A(3,3) // 2 1 0 3 - 
               + A(0,2) * A(1,0) * A(2,1) * A(3,3) // 2 0 1 3 + 
               
               - A(0,0) * A(1,1) * A(2,3) * A(3,2) // 0 1 3 2 - 
               + A(0,0) * A(1,2) * A(2,3) * A(3,1) // 0 2 3 1 + 
               + A(0,1) * A(1,0) * A(2,3) * A(3,2) // 1 0 3 2 + 
               - A(0,1) * A(1,2) * A(2,3) * A(3,0) // 1 2 3 0 - 
               + A(0,2) * A(1,1) * A(2,3) * A(3,0) // 2 1 3 0 + 
               - A(0,2) * A(1,0) * A(2,3) * A(3,1) // 2 0 3 1 - 
               
               + A(0,0) * A(1,3) * A(2,1) * A(3,2) // 0 3 1 2 + 
               - A(0,0) * A(1,3) * A(2,2) * A(3,1) // 0 3 2 1 - 
               - A(0,1) * A(1,3) * A(2,0) * A(3,2) // 1 3 0 2 - 
               + A(0,1) * A(1,3) * A(2,2) * A(3,0) // 1 3 2 0 + 
               - A(0,2) * A(1,3) * A(2,1) * A(3,0) // 2 3 1 0 - 
               + A(0,2) * A(1,3) * A(2,0) * A(3,1) // 2 3 0 1 + 
               
               - A(0,3) * A(1,0) * A(2,1) * A(3,2) // 3 0 1 2 - 
               + A(0,3) * A(1,0) * A(2,2) * A(3,1) // 3 0 2 1 + 
               + A(0,3) * A(1,1) * A(2,0) * A(3,2) // 3 1 0 2 + 
               - A(0,3) * A(1,1) * A(2,2) * A(3,0) // 3 1 2 0 - 
               + A(0,3) * A(1,2) * A(2,1) * A(3,0) // 3 2 1 0 + 
               - A(0,3) * A(1,2) * A(2,0) * A(3,1) // 3 2 0 1 - 
               ;
        
    } else {
      
        Float ret = 0.;
        int sign  = 1;
        
        int i = 77;
        std::vector<int>  aux( A.getdimin() );
        std::vector<int> perm( A.getdimin() );
        for( int j = 0; j < perm.size(); j++ ) perm[j] = j;
        
        HeapsAlgorithmInit( i, aux, perm );
        
        do {
          
          Float summand = sign;
          for( int r = 0; r < A.getdimout(); r++ )
            summand *= A( r, perm[r] );
          
          ret += summand;
            
          sign = -sign;
          
        } while ( HeapsAlgorithmStep( i, aux, perm ) );
        
        return ret;
    }
}


Float Determinant_laplaceexpansion( const DenseMatrix& A )
{
    assert( A.issquare() );
    
    if( A.getdimin() == 0 )
        return 0.;
    
    Float ret = 0.;
    int sign  = 1;
    
    int i = 77;
    std::vector<int>  aux( A.getdimin() );
    std::vector<int> perm( A.getdimin() );
    for( int j = 0; j < perm.size(); j++ ) perm[j] = j;
    
    HeapsAlgorithmInit( i, aux, perm );
    
    do {
      
      Float summand = sign;
      for( int r = 0; r < A.getdimout(); r++ )
        summand *= A( r, perm[r] );
      
      ret += summand;
        
      sign = -sign;
      
    } while ( HeapsAlgorithmStep( i, aux, perm ) );
    
    return ret;
    
}




DenseMatrix CofactorMatrix( const DenseMatrix& A )
{
  assert( A.issquare() );
  
  if( A.getdimin() == 0 ) 
    return DenseMatrix( 0, 0 );
  
  int i = 77;
  std::vector<int>  aux( A.getdimin() - 1 );
  std::vector<int> perm( A.getdimin() - 1 );
  for( int j = 0; j < perm.size(); j++ ) perm[j] = j;
  
  DenseMatrix cof( A.getdimin(), A.getdimin(), 0. );
  
  HeapsAlgorithmInit( i, aux, perm );
  
  int sign_perm = 1;
  
  do {
    
    for( int r = 0; r < A.getdimin(); r++ )
    for( int c = 0; c < A.getdimin(); c++ )
    {
      
      int sign_entry = integerpower( -1, r+c );
      
      Float summand = sign_perm * sign_entry;
      assert( summand == 1. || summand == -1. );
      
      for( int i = 0; i < A.getdimin() - 1; i++ )
        summand *= A( i < c ? i : i+1, perm[i] < r ? perm[i] : perm[i]+1 );
      
      cof(r,c) += summand;
      
    }
    
    sign_perm *= -1;
    
  } while ( HeapsAlgorithmStep( i, aux, perm ) );
  
  return cof;
}


DenseMatrix Inverse( const DenseMatrix& A )
{
    assert( A.issquare() ); 
    return CofactorMatrix( A ) / Determinant( A );
}


void InverseAndDeterminant( const DenseMatrix& A, DenseMatrix& Ainv, Float& Adet )
{
    assert( A.issquare() ); 
    DenseMatrix Cof = CofactorMatrix( A );
    Adet = Determinant( A );
    Ainv = Cof / Adet;
    
}





DenseMatrix SubdeterminantMatrixSquare( const DenseMatrix& A, int k )
{
    A.check();
    assert( A.issquare() );
    assert( 0 <= k && k <= A.getdimin() );
    
    const int n = A.getdimin();
    IndexRange fromrange = IndexRange( 0, k-1 );
    IndexRange torange = IndexRange( 0, n-1 );
    std::vector<IndexMap> sigmas = generateSigmas( fromrange, torange );
    
    DenseMatrix ret( sigmas.size() );
    for( int rim = 0; rim < sigmas.size(); rim++ )
    for( int cim = 0; cim < sigmas.size(); cim++ )
    {
        ret(rim,cim) = determinant( A.submatrix( sigmas.at(rim), sigmas.at(cim) ) );
    }
    
    ret.check();
    return ret;
}


DenseMatrix SubdeterminantMatrix( const DenseMatrix& A, int k )
{
    A.check();
    assert( 0 <= k && k <= A.getdimin() && k <= A.getdimout() );
    
    IndexRange range_from = IndexRange( 0, k-1 );
    IndexRange range_rows = IndexRange( 0, A.getdimout()-1 );
    IndexRange range_cols = IndexRange( 0, A.getdimin()-1 );
    std::vector<IndexMap> sigmas_rows = generateSigmas( range_from, range_rows );
    std::vector<IndexMap> sigmas_cols = generateSigmas( range_from, range_cols );
    
    DenseMatrix ret( sigmas_rows.size(), sigmas_cols.size(), 0. );
    for( int rim = 0; rim < sigmas_rows.size(); rim++ )
    for( int cim = 0; cim < sigmas_cols.size(); cim++ )
    {
        ret(rim,cim) = determinant( A.submatrix( sigmas_rows.at(rim), sigmas_cols.at(cim) ) );
    }
    
    ret.check();
    return ret;
}






