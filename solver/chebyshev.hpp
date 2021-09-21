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
    
    Float* residual_old = residual;
    Float* residual_new = new (std::nothrow) Float[N];
    
    Float* daimonic_old = new (std::nothrow) Float[N];
    Float* daimonic_new = new (std::nothrow) Float[N];
    
    assert( residual_new );
    assert( daimonic_old );
    assert( daimonic_new );
    
    assert( 0. <= lower and lower <= upper );

    const Float gamma = (upper - lower) / 2.0;
    const Float delta = (upper + lower) / 2.0;
    Float alpha = 2. / delta;
    Float beta  = 0.;
    
    Float r_r = 0.;
    
    int K = 0;
    
    while( K < 10 * N ){
        
        bool restart_condition = ( K == 0 ); // or K % 1000 == 0;
        
        bool residual_seems_small = std::sqrt(r_r) < allowed_error;

        if( restart_condition or residual_seems_small ) {
            
            r_r = 0.;
            
            for( int c = 0; c < N; c++ ) {
                
                residual_old[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual_old[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                daimonic_old[c] = 0.;
                
                r_r += residual_old[c] * residual_old[c];
                
            }

            alpha = 2. / delta;

            beta = 0.;
        
        }
        
        /* Check whether residual is small */
                
        bool residual_is_small = std::sqrt(r_r) < allowed_error;
        
        if( residual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        r_r = 0.;

        for( int c = 0; c < N; c++ )
        {
            
            daimonic_new[c] = precon[c] * residual_old[c] + beta * daimonic_old[c];

            x[c] = x[c] + alpha * daimonic_new[c];

            residual_new[c] = residual_old[c];
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                residual_new[c] -= alpha * csrvalues[ d ] * (
                                               precon[csrcolumns[d]] * residual_old[csrcolumns[d]] + beta * daimonic_old[csrcolumns[d]]
                                           );
            
            r_r += residual_new[c] * residual_new[c];
            
        }

        beta = square( gamma * alpha ) / 4.;
        alpha = 1. / ( delta - beta );
        
        std::swap( residual_new, residual_old );
        std::swap( daimonic_new, daimonic_old );

        LOG << K << space << alpha << space << beta << nl;
        
        
        if( print_modulo > 0 and K % print_modulo == 0 ) 
            printf("Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", K, N, (long double)std::sqrt(r_r), (long double) allowed_error );
        
        K++;
        
    }
    
    printf("Final residual after %d of max. %d iterations: %.9Le (%.9Le)\n", K, N, (long double)std::sqrt(r_r), (long double) allowed_error );

    
    if( residual == residual_old ) {
        delete[] residual_new;
    } else { 
        delete[] residual_old;
    }
    delete[] ( daimonic_old );
    delete[] ( daimonic_new );

    

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

    Float r_size = r.norm() / normb;

    if ( r_size < threshold )
        return;

    int recent_iteration_count = 0;
    while( recent_iteration_count < max_iteration_count ) 
    {
        if( r_size < threshold )
            break;                     // convergence
        
        z = invprecon * r;                 // apply preconditioner

        if( recent_iteration_count == 0 ) {
        
            alpha = 2.0 / delta;
            p = z;
        
        } else {
        
            Float beta = square( gamma * alpha / 2.0 );       // calculate new beta
            alpha = 1.0 / ( delta - beta );     // calculate new alpha
            p = z + beta * p;             // update search direction
        
            LOG << recent_iteration_count << space << alpha << space << beta << nl;

        }

        q = A * p;
        x += alpha * p;                 // update approximation vector
        r -= alpha * q;                 // compute residual

        r_size = r.norm() / normb;

        

        recent_iteration_count++;
    }

    LOG << "Final result after " << recent_iteration_count << " of max. " << max_iteration_count << " iterations: " 
            << r_size << "(" << threshold << ")" << nl;
}



#endif
