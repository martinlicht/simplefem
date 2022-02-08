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


inline int CheybyshevIteration_DiagonalPreconditioner( 
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

    LOG << lower << space << upper << nl;
    
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
            LOGPRINTF( "Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", K, N, (long double)std::sqrt(r_r), (long double) allowed_error );
        
        K++;
        
    }
    
    LOGPRINTF( "Final residual after %d of max. %d iterations: %.9Le (%.9Le)\n", K, N, (long double)std::sqrt(r_r), (long double) allowed_error );

    
    delete[] ( x_prev );
    delete[] ( x_curr );
    delete[] ( x_next );

    delete[] ( zaratite );

    return K;

}








#endif
