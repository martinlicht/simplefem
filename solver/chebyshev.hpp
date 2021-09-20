#ifndef INCLUDEGUARD_SOLVER_CHEBYSHEV
#define INCLUDEGUARD_SOLVER_CHEBYSHEV


#include <omp.h>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <utility>

#include "../basic.hpp"






// The Convergence of Inexact Chebyshev and Richardson Iterative Methods for Solving Linear Systems


inline void CheybyshevIteration_DiagonalPreconditioner( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float allowed_error,
    unsigned int print_modulo,
    const Float* __restrict__ precon,
    const Float lower,
    const Float upper
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    assert( allowed_error > 0 );
    assert( print_modulo >= 0 );
    assert( precon );
    
    Float* __restrict__ zaratite = new (std::nothrow) Float[N];
    assert( zaratite );
    
    Float* x_prev = new (std::nothrow) Float[N];
    Float* x_curr = new (std::nothrow) Float[N];
    Float* x_next = new (std::nothrow) Float[N];
    assert( x_prev );
    assert( x_curr );
    assert( x_next );
    
    const Float alpha = 2. / ( upper + lower );
    const Float mu    = ( upper + lower ) / ( upper - lower );
    
    Float gamma_prev = notanumber;
    Float gamma_curr = notanumber;
    Float gamma_next = notanumber;
    
    Float r_r = notanumber;
    
    int K = 0;
    
    while( K < N ){
        
        bool restart_condition = ( K == 0 ); // or K % 1000 == 0;
        
        bool residual_seems_small = std::sqrt(r_r) < allowed_error;

        if( restart_condition or residual_seems_small ) {
            
            gamma_prev = 1.;
            
            gamma_curr = mu;
            
            r_r = 0.;
            
            #if defined(_OPENMP)
            #pragma omp parallel for 
            #endif
            for( int c = 0; c < N; c++ )
                x_prev[c] = x[c];
            
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(+:r_r)
            #endif
            for( int c = 0; c < N; c++ ) {
                
                residual[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] -= csrvalues[ d ] * x_prev[ csrcolumns[d] ];
                
                zaratite[c] = precon[c] * residual[c];
                
                x_curr[c] = x_prev[c] + alpha * zaratite[c];
                
                assert( std::isfinite( x_curr[c] ) );
                
                r_r += residual[c] * residual[c];
                
            }
        
        }
        
        /* Check whether residual is small */
                
        bool residual_is_small = std::sqrt(r_r) < allowed_error;
        
        if( residual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        Float gamma_next = 2. * mu * gamma_curr - gamma_prev; 
        
        Float omega = 2. * mu * gamma_curr / gamma_next; 
        
        r_r = 0.;
            
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:r_r) 
        #endif 
        for( int c = 0; c < N; c++ )
        {
            
            residual[c] = b[c];
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                residual[c] -= csrvalues[ d ] * x_curr[ csrcolumns[d] ];
            
            zaratite[c] = precon[c] * residual[c];
            
            x_next[c] = x_prev[c] + omega * ( alpha * zaratite[c] + x_curr[c] - x_prev[c] );
            
            r_r += residual[c] * residual[c];
            
        }
        
        std::swap( gamma_curr, gamma_prev );
        std::swap( gamma_next, gamma_curr );
        std::swap( x_curr, x_prev );
        std::swap( x_next, x_curr );
        
        
        if( print_modulo > 0 and K % print_modulo == 0 ) 
            printf("Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", K, N, (long double)std::sqrt(r_r), (long double) allowed_error );
        
        K++;
        
    }
    
    printf("Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", K, N, (long double)std::sqrt(r_r), (long double) allowed_error );

    
    delete[] ( x_prev );
    delete[] ( x_curr );
    delete[] ( x_next );

    delete[] ( zaratite );

}






inline void Chebyshev( 
    const LinearOperator& A, 
    const FloatVector& b, 
    FloatVector& x, 
    FloatVector& r, 
    int iterNum, 
    Float lMin, 
    Float lMax 
){
  
  assert( 0. <= lMin and lMin <= lMax );

  const Float delta = (lMax + lMin) / 2.;
  const Float gamma = (lMax - lMin) / 2.;

  r = b - A * x;

  assert( x.isfinite() );
  assert( r.isfinite() );
      
  FloatVector p( x.getdimension() );
  FloatVector z( x.getdimension() );
  Float alpha, beta;
  
  int i;
  for( i = 0; i < iterNum; i++ )
  {
      
      z = r;

      if (i == 0) {

          p = z;
          alpha = 2. / delta;

      } else {
          
          beta = gamma * alpha / 2;
          beta = beta * beta;
          assert( absolute( delta - beta ) > machine_epsilon );
          alpha = 1./( delta - beta );
          p = z + beta * p;

      }

      assert( x.isfinite() );
      assert( r.isfinite() );
      assert( p.isfinite() );
      assert( z.isfinite() );
      

      x = x + alpha * p;
      r = b - A * x; 
      if ( r*r < 1e-15) break; 

      assert( r.isfinite() );

      printf("Cheybyshev Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", i, iterNum, (long double)std::sqrt(r*r), (long double) 1e-15 );

  }

  printf("Cheybyshev Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", i, iterNum, (long double)std::sqrt(r*r), (long double) 1e-15 );

}













//*****************************************************************
// Iterative template routine -- CHEBY
//
// CHEBY solves the symmetric positive definite linear
// system Ax = b using the Preconditioned Chebyshev Method
//
// CHEBY follows the algorithm described on p. 30 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************


void 
CHEBY(const LinearOperator &A, FloatVector &x, const FloatVector &b,
      const LinearOperator& invprecon, 
      const int max_iteration_count, const Float threshold,
      Float eigmin, Float eigmax)
{

    assert( 0. <= eigmin and eigmin <= eigmax );
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    assert( invprecon.getdimin()  == x.getdimension() );
    assert( invprecon.getdimout() == b.getdimension() );
    
    const int dimension = A.getdimin();

    
    Float alpha, beta;
    const Float gamma = (eigmax - eigmin) / 2.0;
    const Float delta = (eigmax + eigmin) / 2.0;
    
    FloatVector p( dimension );
    FloatVector q( dimension );
    FloatVector z( dimension );
    
    FloatVector r = b - A * x;

    Float normb = b.norm();
    if( normb == 0.0 )
        normb = 1;

    Float r_r = r.norm() / normb;

    if ( r_r < threshold*threshold )
        return;

    int recent_iteration_count = 0;
    while( recent_iteration_count < max_iteration_count ) 
    {
        if( r_r < threshold*threshold )
            break;                     // convergence
        
        z = invprecon * r;                 // apply preconditioner

        if( recent_iteration_count == 0 ) {
        
            p = z;
            alpha = 2.0 / delta;
        
        } else {
        
            Float beta = square( gamma * alpha / 2.0 );       // calculate new beta
            alpha = 1.0 / ( delta - beta );     // calculate new alpha
            p = z + beta * p;             // update search direction
        
        }

        q = A * p;
        x += alpha * p;                 // update approximation vector
        r -= alpha * q;                 // compute residual

        r_r = r.norm() / normb;

        

        recent_iteration_count++;
    }

    LOG << "Final result after " << recent_iteration_count << " of max. " << max_iteration_count << " iterations: " 
            << r_r << "(" << threshold*threshold << ")" << nl;
}



#endif
