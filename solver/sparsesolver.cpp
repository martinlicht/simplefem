


#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <memory.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>



#include "sparsesolver.hpp"




void ConjugateGradientSolverCSR( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float allowed_error,
    unsigned int restart_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    assert( allowed_error > 0 );
    assert( restart_modulo > 0 );
    
    Float* __restrict__ direction = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ auxiliary = (Float*)malloc( sizeof(Float) * N );
    assert( direction );
    assert( auxiliary );
    
    Float r_r;

    



    
    // Initialize : compute the residual, the direction, and r * r
    
    r_r = 0.;
    
    #pragma omp parallel for reduction(+:r_r)
    for( int c = 0; c < N; c++ ) {
        
        residual[c] = b[c];
        
        for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
            residual[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
        
        direction[c] = residual[c];
        
        r_r += residual[c] * residual[c];
    }

    

    int K = 0;

    while( K < N ){
        
        
        bool restart_condition = K == 0; // or K % 1000 == 0;
        
        bool residual_seems_small = std::sqrt(r_r) < allowed_error;
        
        if( restart_condition or residual_seems_small ) {
            
            r_r = 0.;
            
            #pragma omp parallel for reduction(+:r_r)
            for( int c = 0; c < N; c++ ) {
                
                residual[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                direction[c] = residual[c]; // this line seems to slow down performance ....
                
                r_r += residual[c] * residual[c];
            }
            
            
        }
        
        /* Check whether residual is small */
        
        bool residual_is_small = std::sqrt(r_r) < allowed_error;
        
        if( residual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        // NOTE The calculation of d_r is reduced to r_r, which is already known.
        
        Float d_r = r_r;
//         Float d_r = 0.;
        Float d_Ad = 0.;
        
        #pragma omp parallel for reduction(+:d_Ad) //d_r,
        for( int c = 0; c < N; c++ )
        {
            //d_r += direction[c] * residual[c];
            
            auxiliary[c] = 0.;
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                auxiliary[c] += csrvalues[ d ] * direction[ csrcolumns[d] ];
                    
            d_Ad += direction[c] * auxiliary[c];
        }
        
        Float alpha = d_r / d_Ad;
        
        Float r_r_new = 0.;
        
        #pragma omp parallel for reduction(+:r_r_new) //r_r_old
        for( int c = 0; c < N; c++ )
        {
            
            x[c] += alpha * direction[c];
            
            residual[c] -= alpha * auxiliary[c];
            
            r_r_new += residual[c] * residual[c];
        }
        
        Float beta = r_r_new / r_r;
        
        r_r = r_r_new;
        
        #pragma omp parallel for
        for( int c = 0; c < N; c++ )
            direction[c] = residual[c] + beta * direction[c];
//             auxiliary[c] = residual[c] + beta * direction[c]; // NOTE: Disabled swapping below, simply forget auxiliary buffer
        
        
        // main part of iteration done. Swap buffers and give output
        
//         {
//             Float* swappi = auxiliary;
//             auxiliary = direction;
//             direction = swappi;
//         }
        
//         if( K % 100 == 0 ) 
//         printf("At Iteration %d we have %.9Le --- [%.9Le,%.9Le,%.9Le,%.9Le]\n",
//             K,
//             (long double)std::sqrt(r_r), (long double)alpha, (long double)beta,
//             (long double)d_Ad, (long double)r_r_new );
        
        K++;
        
        
    };
    
    printf("Residual after %d of max. %d iterations: %.9Le\n", K, N, (long double)std::sqrt(r_r) );

    
    free( direction ); 
    free( auxiliary );

}































// s-step iterative methods for symmetric linear systems. Chronopoulos & Gear



void ConjugateGradientSolverCSR_variant( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ r,
    const Float allowed_error,
    unsigned int restart_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( r );
    assert( allowed_error > 0 );
    assert( restart_modulo > 0 );
    
    Float* __restrict__  d = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ Ad = (Float*)malloc( sizeof(Float) * N );
    assert(  d );
    assert( Ad );
    
    Float alpha;
    Float  beta;
    Float   r_r;

    
    // Initialize : compute the residual, the direction, and r * r
    
    int K = 0;

    while( K < N ){
        
        
        bool restart_condition = (K == 0); // or (K % 1000 == 1);
        
        bool residual_seems_small = (std::sqrt(r_r) < allowed_error);
        
        if( restart_condition or residual_seems_small ) {
            
            r_r = 0.;

            #pragma omp parallel for reduction(+:r_r)
            for( int c = 0; c < N; c++ ) {
                
                r[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    r[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                d[c] = r[c]; 
                
                r_r += r[c] * r[c];
                
            }
            
            Float Ad_d = 0.;
            
            #pragma omp parallel for reduction(+:Ad_d)
            for( int c = 0; c < N; c++ ) {
                
                Ad[c] = 0.;
                
                for( int e = csrrows[c]; e < csrrows[c+1]; e++ )
                    Ad[c] += csrvalues[ e ] * d[ csrcolumns[e] ];
                
                Ad_d += Ad[c] * d[c];
                
            }
            
            alpha = r_r / Ad_d;
            
            beta = 0.;
            
        }
        
        /* Check whether residual is small */
        
        bool residual_is_small = std::sqrt(r_r) < allowed_error;
        
        if( residual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        
        Float r_r_new = 0.;
        
        #pragma omp parallel for reduction(+:r_r_new)
        for( int c = 0; c < N; c++ )
        {
            
            // update x and r with the current direction
            
            x[c] = x[c] + alpha *  d[c];
            
            r[c] = r[c] - alpha * Ad[c];
            
            r_r_new += r[c] * r[c];
            
            // update the direction 
            
            d[c] = r[c] + beta * d[c];
            
            
        }
        
        Float Ad_d  = 0.;
        Float Ad_Ad = 0.;
        
        #pragma omp parallel for reduction(+:Ad_d,Ad_Ad) //r_r_old
        for( int c = 0; c < N; c++ )
        {
            
            Ad[c] = 0.;
                
            for( int e = csrrows[c]; e < csrrows[c+1]; e++ )
                Ad[c] += csrvalues[ e ] * d[ csrcolumns[e] ];
            
            Ad_d  += Ad[c] *  d[c];
            
            Ad_Ad += Ad[c] * Ad[c];
        }
        
//         assert( Ad_d    >= 0 );
//         assert( Ad_Ad   >= 0 );
//         assert( r_r_new >= 0 );
        
        
        alpha = r_r_new / Ad_d;
        
        beta = ( alpha * alpha * Ad_Ad - r_r_new ) / ( r_r_new );
        
        r_r = r_r_new;
                
        
//         if( K % 100 == 0 ) 
//         printf("At Iteration %d we have %.9Le --- [%.9Le,%.9Le,%.9Le,%.9Le]\n",
//             K,
//             (long double)std::sqrt(r_r), (long double)alpha, (long double)beta,
//             (long double)Ad_d, (long double)Ad_Ad );
        
        K++;
        
        
    };
    
    printf("Residual after %d of max. %d iterations: %.9Le\n", K, N, (long double)std::sqrt(r_r) );

    
    free(  d ); 
    free( Ad );

}














































// https://pdf.sciencedirectassets.com/271610/1-s2.0-S0377042700X01380/1-s2.0-0377042789900459/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEEYaCXVzLWVhc3QtMSJHMEUCIQC0fGW%2BgRGDBWFchmZEoI9WEfGcrraFFHFlxAvZcabrCwIgGRxFWZ0FrfCLC4jOxtmY4OqCj%2BM88HmkPN5HvYkN0ywqtAMIHxADGgwwNTkwMDM1NDY4NjUiDEic7xJsav6PnYjOViqRA1pP5BCisC3iJlcLEsGlD7nomRKLfUm11vEYQJBRoPi%2Fus1pCt0Imtxp7i6KUiPyhqHLCYlvwkUHHjjsUI0cNiKnay88cVt%2F9ASRgRAEfv70P0hOpm9M7XCA%2BwOIbua%2B%2B4JwRzjyBRc4FX9n5ubOuRhya3OsB%2BCw8Bzl3AdbsMeILAAktCfoS%2FLjO94yJvxg8%2FQZ4I%2BCIujpqDIGAPyNVawq7eyf4gtW8hAL12BQnVKRw6N9IDK1eHHkO5bUbA4Otzc5YhzUYkOGNx08B6Y%2BaXe33fkjRIlAwZJTDR0bD%2FN88eCYDwnjDTLEv%2FRlgJRxO5X%2BQJ6su0%2F8gUc6QBW3HJV6WRdlJxIl98nFKRQcSrQ9%2BO9lANfoshHF24VF%2FbszqHA0NDqDCovATdWrl557mjeEBTjthYiICKKlF7pDGDDr4rJgR6i9nY0toDvzm%2BoWV%2Fx8a08W3Ps3FOhD4GO1MQ9egc6OVItC%2BPjDgO5D4iJMx2%2F6hIOZw1rbM8tfj2FUuJmTHBosXTTFUkqlud%2F%2FRvCSMN23vPkFOusBLJH1dKnMDj5dtKJIqK8ldhkUiP7rmuwXaXVfuCukRLbM2%2F6sgR%2BBRZt04K08oWUomBV0tqSmG9XN%2BwVphEey1sqgnKg%2Fs9E7N8z6Dkw5feeplyssnN%2FTnBYlN7aoOTPrkh1PRaWQnku%2FK77fT9T0ImU7FWwG5yszdtqkG19UCUjE%2BEwQ94XHoQ8BGHs0roSKrnl%2ByminvDBr1KW7hMToKJ2lNIx6eyL%2BPzQe45aaSE1cD%2FAgdK3jbmIJ1Ba%2BaMIabqf2QcJsspWh2Pfah9cU2EGdSdcC4jeOVvzpWi1Vj%2Bk3DVUjAUxeejkqXA%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20200808T233708Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYUWZMLHGF%2F20200808%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=d2bbeca441d5735625e00f26ea4a38be9bd96a1a0e933e51611bf07d0add0df1&hash=40f879e2028f81408eb3c5365387249adb2c84662044f49bc6aed628fe0ea4f8&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=0377042789900459&tid=spdf-a4149c9d-4592-44d8-9c0c-8acec3388eb3&sid=0182abe95a749047995ba14457d108682c17gxrqb&type=client

void ConjugateGradientSolverCSR_sstep( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ r,
    const Float allowed_error,
    unsigned int restart_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( r );
    assert( allowed_error > 0 );
    assert( restart_modulo > 0 );
    
    Float* __restrict__ d = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ Ad = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ Ar = (Float*)malloc( sizeof(Float) * N );
    assert(  d );
    assert( Ad );
    assert( Ar );
    
    Float alpha;
    Float  beta;
    Float   r_r;

    
    // Initialize : compute the residual, the direction, and r * r
    
    int K = 0;

    while( K < N ){
        
        
        bool restart_condition = (K == 0); // or (K % 1000 == 1);
        
        bool residual_seems_small = (std::sqrt(r_r) < allowed_error);
        
        if( restart_condition or residual_seems_small ) {
            
            r_r = 0.;

            #pragma omp parallel for reduction(+:r_r)
            for( int c = 0; c < N; c++ ) {
                
                r[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    r[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                d[c] = r[c]; 
                
                r_r += r[c] * r[c];
                
            }
            
            Float d_r = r_r;

            Float Ad_d = 0.;
            
            #pragma omp parallel for reduction(+:Ad_d)
            for( int c = 0; c < N; c++ ) {
                
                Ar[c] = 0.;
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    Ar[c] += csrvalues[ d ] * r[ csrcolumns[d] ];
                
                Ad[c] = Ar[c];
                
                Ad_d += Ad[c] * d[c];
                
            }
            
            alpha = d_r / Ad_d;
            
            beta = 0.;
            
        }
        
        /* Check whether residual is small */
        
        bool residual_is_small = std::sqrt(r_r) < allowed_error;
        
        if( residual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        
        Float r_r_new = 0.;
        
//         Float Ad_d = 0.;
//         Float d_r = 0.;
        
        #pragma omp parallel for reduction(+:r_r_new)
        for( int c = 0; c < N; c++ )
        {
            
            // compute the step
            
            d[c] = r[c] + beta * d[c];
            
            Ad[c] = Ar[c] + beta * Ad[c];
            
            // go along that step 
            
            x[c] = x[c] + alpha * d[c];
            
            r[c] = r[c] - alpha * Ad[c];
            
            r_r_new += r[c] * r[c];
            
//             Ad_d += Ad[c] * d[c];
//             
//             d_r += d[c] * r[c];
        }
        
        Float Ar_r = 0.;
        
        #pragma omp parallel for reduction(+:Ar_r) //r_r_old
        for( int c = 0; c < N; c++ )
        {
            
            Ar[c] = 0.;
                
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                Ar[c] += csrvalues[ d ] * r[ csrcolumns[d] ];
            
            Ar_r += r[c] * Ar[c];
            
        }
        
//         assert( r_r > 0 );
//         assert( Ad_d > 0 );
        
        
        beta = r_r_new / r_r;
        
        alpha = r_r_new / ( Ar_r - r_r_new * beta / alpha );
        
        r_r = r_r_new;
                
        
//         if( K % 100 == 0 ) 
//         printf("At Iteration %d we have %.9Le --- [%.9Le,%.9Le,%.9Le,%.9Le]\n",
//             K,
//             (long double)std::sqrt(r_r), (long double)alpha, (long double)beta,
//             (long double)Ad_d, (long double)d_r );
        
        K++;
        
        
    };
    
    printf("Residual after %d of max. %d iterations: %.9Le\n", K, N, (long double)std::sqrt(r_r) );

    
    free(  d ); 
    free( Ad );
    free( Ar );

}































void ConjugateResidualSolverCSR( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ r,
    const Float allowed_error,
    unsigned int restart_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( r );
    assert( allowed_error > 0 );
    assert( restart_modulo > 0 );
    
    Float* __restrict__  d = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ Ad = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ Ar = (Float*)malloc( sizeof(Float) * N );
    assert(  d );
    assert( Ad );
    assert( Ar );
    
    Float alpha;
    Float  beta;
    Float  Ar_r = 0.;

    
    // Initialize : compute the residual, the direction, and r * r
    
    int K = 0;

    while( K < N ){
        
        
        bool restart_condition = (K == 0); // or (K % 1000 == 1);
        
        bool residual_seems_small = (std::sqrt(Ar_r) < allowed_error);
        
        if( restart_condition or residual_seems_small ) {
            
            #pragma omp parallel for
            for( int c = 0; c < N; c++ ) {
                
                r[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    r[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                d[c] = r[c]; 
                
            }
            
            Ar_r = 0.;

            Float Ad_Ad = 0.;
            
            #pragma omp parallel for reduction(+:Ar_r,Ad_Ad)
            for( int c = 0; c < N; c++ ) {
                
                Ar[c] = 0.;
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    Ar[c] += csrvalues[ d ] * r[ csrcolumns[d] ];
                
                Ad[c] = Ar[c];
                
                Ar_r += Ar[c] * r[c];
                
                Ad_Ad += Ad[c] * Ad[c];
                
            }
            
            Float Ad_r = Ar_r;

            alpha = Ad_r / Ad_Ad;
            
            beta = 0.;
            
        }
        
        /* Check whether residual is small */
        
        bool residual_is_small = std::sqrt(Ar_r) < allowed_error;
        
        if( residual_is_small )
            break;

//         std::cout << N << space << Ar_r << space << std::sqrt(Ar_r) << nl;
        

        /* now the main work of the entire algorithm */
        
        
//         Float Ad_Ad = 0.;
//         Float Ad_r = 0.;
        
        #pragma omp parallel for
        for( int c = 0; c < N; c++ )
        {
            
            // compute the step
            
            d[c] = r[c] + beta * d[c];
            
            Ad[c] = Ar[c] + beta * Ad[c];
            
            // go along that step 
            
            x[c] = x[c] + alpha * d[c];
            
            r[c] = r[c] - alpha * Ad[c];
            
//             Ad_Ad += Ad[c] * Ad[c];
//             
//             Ad_r += Ad[c] * r[c];
            
        }
        
//         std::cout << Ar_r << nl;
        
//         assert( debug_r_r > 0 );
        
        Float Ar_r_new = 0.;
        
        Float Ar_Ar = 0.;
        
        #pragma omp parallel for reduction(+:Ar_r_new,Ar_Ar) //r_r_old
        for( int c = 0; c < N; c++ )
        {
            
            Ar[c] = 0.;
                
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                Ar[c] += csrvalues[ d ] * r[ csrcolumns[d] ];
        
            Ar_r_new += Ar[c] *  r[c];
            
            Ar_Ar    += Ar[c] * Ar[c];
            
        }
        
        
        
        beta = Ar_r_new / Ar_r;
        
//         alpha = Ar_r_new / Ad_Ad; // 
        alpha = Ar_r_new / ( Ar_Ar - beta / alpha * Ar_r_new );
        
        
//         if( K % 100 == 0 ) 
//         printf("At Iteration %d we have %.9Le --- [%.9Le,%.9Le,%.9Le,%.9Le]\n",
//             K,
//             (long double)std::sqrt(Ar_r), (long double)alpha, (long double)beta,
//             (long double)Ad_r, (long double)Ar_r_new );
        
//         assert( Ar_r     > 0 );
//         assert( Ar_r_new > 0 );
//         assert( Ad_Ad    > 0 );

        Ar_r = Ar_r_new;
                
        K++;
        
        
    };
    
    printf("Residual after %d of max. %d iterations: %.9Le\n", K, N, (long double)std::sqrt(Ar_r) );

    
    free(  d ); 
    free( Ad );
    free( Ar );

}









void ConjugateResidualSolverCSR_old( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float allowed_error,
    unsigned int restart_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    assert( allowed_error > 0 );
    
    Float* __restrict__ direction = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ auxiliary = (Float*)malloc( sizeof(Float) * N );
    assert( direction );
    assert( auxiliary );
    
    Float Ar_r = 0.;
    
//     #pragma omp parallel for reduction(+:Ar_r)
//     for( int c = 0; c < N; c++ ) {
//         
//         residual[c] = b[c];
//         
//         for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
//             residual[c] += -1. * csrvalues[ d ] * x[ csrcolumns[d] ];
//         
//         direction[c] = residual[c];
//         
//         Ar_r += residual[c] * residual[c];
//         
//     }
    
    int K = 0;
    
    while( K < N ){
                
        bool restart_condition = (K == 0); // or (K % 1000 == 1);
        
        bool residual_seems_small = (std::sqrt(Ar_r) < allowed_error);
        
        if( restart_condition or residual_seems_small ) {
        
            #pragma omp parallel for reduction(+:Ar_r)
            for( int c = 0; c < N; c++ ) {
                
                residual[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] += -1. * csrvalues[ d ] * x[ csrcolumns[d] ];
                
                direction[c] = residual[c];
                
            }
            
            Ar_r = 0.;
            
            #pragma omp parallel for reduction(+:Ar_r)
            for( int c = 0; c < N; c++ ) {
                
                Float Ar_c = 0.;
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    Ar_c += csrvalues[ d ] * residual[ csrcolumns[d] ];
                
                Ar_r += residual[c] * Ar_c;
                
            }
            
        }
        
        
        
        /* Check whether residual is small */
        bool energy_is_small = std::sqrt(Ar_r) < allowed_error;
        
        if( energy_is_small )
            break;

        
        Float r_Ad = 0.;
        
        Float Ad_Ad = 0.;
        
        for( int c = 0; c < N; c++ ) {
            
            Float Ad_c = 0.;
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                Ad_c += csrvalues[ d ] * direction[ csrcolumns[d] ];
            
            r_Ad += residual[c] * Ad_c;
            
            Ad_Ad += Ad_c * Ad_c;
        
        }
        
        Float alpha = r_Ad / Ad_Ad;
        
        for( int c = 0; c < N; c++ ) {
            
            x[c] += alpha * direction[c];
        
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                residual[c] -= alpha * csrvalues[ d ] * direction[ csrcolumns[d] ];
            
        }
        
        Float Ar_r_new = 0.;
        for( int c = 0; c < N; c++ ) {
            
            Float Ar_c = 0.;
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                Ar_c += csrvalues[ d ] * residual[ csrcolumns[d] ];
            
            Ar_r_new += residual[c] * Ar_c;
            
        }
        
        Float beta = Ar_r_new / r_Ad;
        
        for( int c = 0; c < N; c++ ) {
            
            direction[c] = residual[c] + beta * direction[c];
        }
        
        
//         printf("At Iteration %d we have %.9Le --- [%.9Le,%.9Le,%.9Le,%.9Le]\n", K,
//                     (long double)sqrt(Ar_r_new), (long double)alpha, (long double)beta,
//                     (long double)r_Ad, (long double)Ar_r );
        
//         if( std::sqrt( Ar_r_new) < allowed_error ) 
//             break;
        
        Ar_r = Ar_r_new;
        
        K++;
    };
    
    printf("Residual after %d of max. %d iterations: %.9Le\n", K, N, (long double)sqrt(Ar_r) );
    

    free( direction ); 
    free( auxiliary );

}








                      





// void UzawaConjugateResidualCSR( const int M, const int pdim,
//                                 Float* const sigma, Float* const u,
//                                 Float* const ressigma, Float *const resu,
//                                 const Float* const e, const Float* const f,
//                                 int nShadowedSigma, int nShadowedU,
//                                 const int* shadowSigma, const int* shadowU,
//                                 const Float* const entriesA, const int* const csrrowsA, const int* const csrcolumnsA,
//                                 const Float* const entriesB, const int* const csrrowsB, const int* const csrcolumnsB,
//                                 const Float* const entriesC, const int* const csrrowsC, const int* const csrcolumnsC
//                             )
// {
//     
//     /**************************
//       A Bt | s  =  e
//       B C  | u  =  f
//     ***************************
//       exactly this sign
//       convention!!!
//     **************************/
//     assert( M > 0 );
//     assert( pdim > 0 );
//     assert( sigma );
//     assert( u );
//     assert( e );
//     assert( f );
//     assert( nShadowedSigma >= 0 );
//     assert( nShadowedU >= 0 );
//     if( nShadowedSigma > 0 ) assert( shadowSigma );
//     if( nShadowedU > 0 ) assert( shadowU );
//     assert( entriesA ); assert( csrrowsA ); assert( csrcolumnsA );
//     assert( entriesB ); assert( csrrowsB ); assert( csrcolumnsB );
//     assert( entriesC ); assert( csrrowsC ); assert( csrcolumnsC );
//     
//     Float* entriesBt = NULL;
//     int*   csrrowsBt = NULL;
//     int*   csrcolumnsBt = NULL;
//     
//     sparsetranspose( pdim, M,
//        (const int*)csrrowsB, (const int*)csrcolumnsB, (const Float*)entriesB,
//        &csrrowsBt, &csrcolumnsBt, &entriesBt );
//     
// //     printsparse( M, csrrowsBt, csrcolumnsBt, entriesBt );
//     
//     check_matrix( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt );
//     check_matrix( pdim, M, csrrowsB,  csrcolumnsB,  entriesB );
//     check_matrix( pdim, pdim, csrrowsC, csrcolumnsC, entriesC );
//     check_matrix( M, M, csrrowsA, csrcolumnsA, entriesA );
//     
//     Float* tempM = (Float*)mallocfloat(    M ); assert(    M );
//     Float* tempN = (Float*)mallocfloat( pdim ); assert( pdim );
//     
//     setfloats( M,    ressigma, 0. );
//     setfloats( pdim, resu,     0. );
//     setfloats( M,    tempM,    0. );
//     setfloats( pdim, tempN,    0. );
//     
//     Float* r       = (Float*)mallocfloat( pdim ); assert(       r );
//     Float* C_r     = (Float*)mallocfloat( pdim ); assert(     C_r );
//     Float* AiBt_r  = (Float*)mallocfloat(    M ); assert(  AiBt_r );
//     Float* BAiBt_r = (Float*)mallocfloat( pdim ); assert( BAiBt_r );
//     Float* d       = (Float*)mallocfloat( pdim ); assert(       d );
//     Float* C_d     = (Float*)mallocfloat( pdim ); assert(     C_d );
//     Float* AiBt_d  = (Float*)mallocfloat(    M ); assert(  AiBt_d );
//     Float* BAiBt_d = (Float*)mallocfloat( pdim ); assert( BAiBt_d );
//     Float* p       = (Float*)mallocfloat( pdim ); assert(       p );
//     Float* C_p     = (Float*)mallocfloat( pdim ); assert(     C_p );
//     Float* AiBt_p  = (Float*)mallocfloat(    M ); assert(  AiBt_p );
//     Float* BAiBt_p = (Float*)mallocfloat( pdim ); assert( BAiBt_p );
// 
//     setfloats( pdim,       r, 0. );
//     setfloats( pdim,     C_r, 0. );
//     setfloats(    M,  AiBt_r, 0. );
//     setfloats( pdim, BAiBt_r, 0. );
//     setfloats( pdim,       d, 0. );
//     setfloats( pdim,     C_d, 0. );
//     setfloats(    M,  AiBt_d, 0. );
//     setfloats( pdim, BAiBt_d, 0. );
//     setfloats( pdim,       p, 0. );
//     setfloats( pdim,     C_p, 0. );
//     setfloats(    M,  AiBt_p, 0. );
//     setfloats( pdim, BAiBt_p, 0. );
//     
//     int maxiter = 10 * (M + pdim);
//     int iter = 0;
//     
//     const Float allowed_error = allowed_error;
//     
//     /*************/
//     /* MAIN LOOP */
//     /*************/
//     
//     while( iter < maxiter ) {
//         
//         printf("@@@@@@@@@@ Uzawa-CRM Iteration %d / %d: unorm is %f\n", iter, maxiter, vectornorm(pdim,u) );
//         
//         /* Start or Restart condition check */
//         if( iter == 0 or ( false && iter % 1000 == 0 ) ) {
//            
//             // tempM = B^t u, tempM = -tempM, tempM = tempM + e
//             sparsematrixvectormultiply( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                                         (const Float*)u, tempM, nShadowedU, shadowU );
//             
//             FREEZE( tempM, nShadowedSigma, shadowSigma );
//                                         
//             scale( tempM, M, -1. );
//             
//             add( tempM, (const Float*)e, M );
//             
//             FREEZE( tempM, nShadowedSigma, shadowSigma );
//             
//             // sigma = A^\inv tempM
//             setfloats( M, sigma, 0. );
//             ConjugateResidualSolverCSR( M, sigma, tempM, csrrowsA, csrcolumnsA, entriesA,
//                                         ressigma, nShadowedSigma, shadowSigma );
//             
//             FREEZE( sigma, nShadowedSigma, shadowSigma );
//             
//             // r = d = f - B sigma - C u
//             copyfloats( r, f, pdim );
//             
//             sparsematrixvectormultiply( pdim, M, csrrowsB, csrcolumnsB, entriesB,
//                                         (const Float*)sigma, tempN, nShadowedSigma, shadowSigma );
//                                         
//             addscaled( r, (const Float*)tempN, pdim, -1. );
//             
//             sparsematrixvectormultiply( pdim, pdim, csrrowsC, csrcolumnsC, entriesC,
//                                         (const Float*)u, tempN, nShadowedU, shadowU );
//                                         
//             addscaled( r, (const Float*)tempN, pdim, -1. );
//             
//             FREEZE( r, nShadowedU, shadowU );
//             
//             // residual Auxiliary variables
//             sparsematrixvectormultiply( pdim, pdim,  csrrowsC,  csrcolumnsC,  entriesC,
//                                              (const Float*)r,     C_r,     nShadowedU,     shadowU );
//                                              
//             FREEZE( C_r, nShadowedU, shadowU );
//             
//             sparsematrixvectormultiply(    M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                                              (const Float*)r,   tempM,     nShadowedU,     shadowU );
//             
//             FREEZE( tempM, nShadowedSigma, shadowSigma );
//             
//             setfloats( M, AiBt_r, 0. );
//             ConjugateResidualSolverCSR( M, AiBt_r, tempM, csrrowsA, csrcolumnsA, entriesA,
//                                         ressigma, nShadowedSigma, shadowSigma );
//                                         
//             FREEZE( AiBt_r, nShadowedSigma, shadowSigma );
//             
//             sparsematrixvectormultiply( pdim,    M,  csrrowsB,  csrcolumnsB,  entriesB,
//                                         (const Float*)AiBt_r, BAiBt_r, nShadowedSigma, shadowSigma );
//             
//             FREEZE( BAiBt_r, nShadowedU, shadowU );
//             
//             // direction auxiliary variables
//             cpyfloats( pdim,       d, (const Float*)r       );
//             cpyfloats( pdim,     C_d, (const Float*)C_r     );
//             cpyfloats(    M,  AiBt_d, (const Float*)AiBt_r  );
//             cpyfloats( pdim, BAiBt_d, (const Float*)BAiBt_r );
//             
//         }
// 
//         Float residualnorm = vectornorm( pdim, r );
// 
//         printf( "@@@@@@@@@@ Residual norm: %f\n", residualnorm );
//         
//         if( residualnorm < allowed_error ) {
//             printf("@@@@@@@@@@ Threshold deceeded.\n");
//             break;
//         }
//         
//         FREEZE(        u,     nShadowedU,     shadowU );
//         FREEZE(    sigma, nShadowedSigma, shadowSigma );
//         FREEZE(     resu,     nShadowedU,     shadowU );
//         FREEZE( ressigma, nShadowedSigma, shadowSigma );
//         FREEZE(    tempM, nShadowedSigma, shadowSigma );
//         FREEZE(    tempN,     nShadowedU,     shadowU );
//         FREEZE(        r,     nShadowedU,     shadowU );
//         FREEZE(      C_r,     nShadowedU,     shadowU );
//         FREEZE(   AiBt_r, nShadowedSigma, shadowSigma );
//         FREEZE(  BAiBt_r,     nShadowedU,     shadowU );
//         FREEZE(        d,     nShadowedU,     shadowU );
//         FREEZE(      C_d,     nShadowedU,     shadowU );
//         FREEZE(   AiBt_d, nShadowedSigma, shadowSigma );
//         FREEZE(  BAiBt_d,     nShadowedU,     shadowU );
//         FREEZE(        p,     nShadowedU,     shadowU );
//         FREEZE(      C_p,     nShadowedU,     shadowU );
//         FREEZE(   AiBt_p, nShadowedSigma, shadowSigma );
//         FREEZE(  BAiBt_p,     nShadowedU,     shadowU );
//         
//         
//         
//         /* compute first variables */
//         
//         // t1
//         Float t1 = scalarproduct( pdim, r, C_r ) - scalarproduct( pdim, r, BAiBt_r );
//         
//         // p, C_p, AiBt_p, BAiBt_p == M d, C M d, AiBt M d, B Ai Bt M d
//         copyfloats( p, (const Float*)C_d,     pdim );
//         sub       ( p, (const Float*)BAiBt_d, pdim );
//         // ...
//         sparsematrixvectormultiply( pdim, pdim,  csrrowsC,  csrcolumnsC,  entriesC,
//                                     (const Float*)p, C_p, nShadowedU, shadowU );
//                                     
//         sparsematrixvectormultiply( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                                     (const Float*)p, tempM, nShadowedU, shadowU );
//         
//         FREEZE(        u,     nShadowedU,     shadowU );
//         FREEZE(    sigma, nShadowedSigma, shadowSigma );
//         FREEZE(     resu,     nShadowedU,     shadowU );
//         FREEZE( ressigma, nShadowedSigma, shadowSigma );
//         FREEZE(    tempM, nShadowedSigma, shadowSigma );
//         FREEZE(    tempN,     nShadowedU,     shadowU );
//         FREEZE(        r,     nShadowedU,     shadowU );
//         FREEZE(      C_r,     nShadowedU,     shadowU );
//         FREEZE(   AiBt_r, nShadowedSigma, shadowSigma );
//         FREEZE(  BAiBt_r,     nShadowedU,     shadowU );
//         FREEZE(        d,     nShadowedU,     shadowU );
//         FREEZE(      C_d,     nShadowedU,     shadowU );
//         FREEZE(   AiBt_d, nShadowedSigma, shadowSigma );
//         FREEZE(  BAiBt_d,     nShadowedU,     shadowU );
//         FREEZE(        p,     nShadowedU,     shadowU );
//         FREEZE(      C_p,     nShadowedU,     shadowU );
//         FREEZE(   AiBt_p, nShadowedSigma, shadowSigma );
//         FREEZE(  BAiBt_p,     nShadowedU,     shadowU );
//         
//         
//         setfloats( M, AiBt_p, 0. );
//         ConjugateResidualSolverCSR( M, AiBt_p, (const Float*)tempM, csrrowsA, csrcolumnsA, entriesA,
//                                     ressigma, nShadowedSigma, shadowSigma );
//                                     
//         sparsematrixvectormultiply( pdim, M, csrrowsB, csrcolumnsB, entriesB,
//                                     (const Float*)AiBt_p, BAiBt_p, nShadowedSigma, shadowSigma );
//         
//         FREEZE(        u,     nShadowedU,     shadowU );
//         FREEZE(    sigma, nShadowedSigma, shadowSigma );
//         FREEZE(     resu,     nShadowedU,     shadowU );
//         FREEZE( ressigma, nShadowedSigma, shadowSigma );
//         FREEZE(    tempM, nShadowedSigma, shadowSigma );
//         FREEZE(    tempN,     nShadowedU,     shadowU );
//         FREEZE(        r,     nShadowedU,     shadowU );
//         FREEZE(      C_r,     nShadowedU,     shadowU );
//         FREEZE(   AiBt_r, nShadowedSigma, shadowSigma );
//         FREEZE(  BAiBt_r,     nShadowedU,     shadowU );
//         FREEZE(        d,     nShadowedU,     shadowU );
//         FREEZE(      C_d,     nShadowedU,     shadowU );
//         FREEZE(   AiBt_d, nShadowedSigma, shadowSigma );
//         FREEZE(  BAiBt_d,     nShadowedU,     shadowU );
//         FREEZE(        p,     nShadowedU,     shadowU );
//         FREEZE(      C_p,     nShadowedU,     shadowU );
//         FREEZE(   AiBt_p, nShadowedSigma, shadowSigma );
//         FREEZE(  BAiBt_p,     nShadowedU,     shadowU );
//         
//         // t2
//         Float t2 = scalarproduct( pdim, d, C_p ) - scalarproduct( pdim, d, BAiBt_p );
//         
//         // alpha
//         Float alpha = t1 / t2;
//         
//         // Update u, sigma
//         addscaled(       u,  (const Float*)     d, pdim,  alpha );
//         addscaled(   sigma,  (const Float*)AiBt_d,    M, -alpha );
//         
//         // Update residual terms
//         addscaled(       r, (const Float*)      p, pdim, -alpha );
//         addscaled(     C_r, (const Float*)    C_p, pdim, -alpha );
//         addscaled(  AiBt_r, (const Float*) AiBt_p,    M, -alpha );
//         addscaled( BAiBt_r, (const Float*)BAiBt_p, pdim, -alpha );
//         
//         // Compute beta
//         Float beta = scalarproduct( pdim, r, C_r ) - scalarproduct( pdim, r, BAiBt_r ); 
//         beta = beta / t1;
//         
//         // Update direction terms
//         scalefloats( pdim,       d, beta );
//         scalefloats( pdim,     C_d, beta );
//         scalefloats(    M,  AiBt_d, beta );
//         scalefloats( pdim, BAiBt_d, beta );
//         
//         add(       d, (const Float*)      r, pdim );
//         add(     C_d, (const Float*)    C_r, pdim );
//         add(  AiBt_d, (const Float*) AiBt_r,    M );
//         add( BAiBt_d, (const Float*)BAiBt_r, pdim );
//         
//         iter++;
//         
//     }
//     
//     /******************/
//     /*** END LOOP *****/
//     /******************/
//     
//     
//     FREEZE(        u,     nShadowedU,     shadowU );
//     FREEZE(    sigma, nShadowedSigma, shadowSigma );
//     FREEZE(     resu,     nShadowedU,     shadowU );
//     FREEZE( ressigma, nShadowedSigma, shadowSigma );
//     FREEZE(    tempM, nShadowedSigma, shadowSigma );
//     FREEZE(    tempN,     nShadowedU,     shadowU );
//     FREEZE(        r,     nShadowedU,     shadowU );
//     FREEZE(      C_r,     nShadowedU,     shadowU );
//     FREEZE(   AiBt_r, nShadowedSigma, shadowSigma );
//     FREEZE(  BAiBt_r,     nShadowedU,     shadowU );
//     FREEZE(        d,     nShadowedU,     shadowU );
//     FREEZE(      C_d,     nShadowedU,     shadowU );
//     FREEZE(   AiBt_d, nShadowedSigma, shadowSigma );
//     FREEZE(  BAiBt_d,     nShadowedU,     shadowU );
//     FREEZE(        p,     nShadowedU,     shadowU );
//     FREEZE(      C_p,     nShadowedU,     shadowU );
//     FREEZE(   AiBt_p, nShadowedSigma, shadowSigma );
//     FREEZE(  BAiBt_p,     nShadowedU,     shadowU );
//     
//     
//     {
//             
//         copyfloats( ressigma, (const Float*)e, M );
//         FREEZE( ressigma, nShadowedSigma, shadowSigma );
//     
//         FREEZE( tempM, nShadowedSigma, shadowSigma );
//         sparsematrixvectormultiply( M, M, csrrowsA, csrcolumnsA, entriesA,
//                                     (const Float*)sigma, tempM, nShadowedSigma, shadowSigma );
//                                     
//         FREEZE( sigma, nShadowedSigma, shadowSigma );
//         sub( ressigma, (const Float*)tempM, M );
//         
//         FREEZE( u, nShadowedU, shadowU );
//         sparsematrixvectormultiply( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                                     (const Float*)u, tempM, nShadowedU, shadowU );
//                                     
//         sub( ressigma, (const Float*)tempM, M );
//         FREEZE( ressigma, nShadowedSigma, shadowSigma );
//         
//         copyfloats( resu, r, pdim );
//         FREEZE( resu, nShadowedU, shadowU );
//         
//         
//     }
//     
// //     printf( "@@@@@@@@@@ Residual sigma %f \n", vectornorm( M, ressigma ) );
// //     printf( "@@@@@@@@@@ Residual u     %f \n", vectornorm( pdim, resu  ) );
//     
//     if(true){
//     
//         // sigma and ressigma
//         // tempM = B^t u, tempM = -tempM, tempM = tempM + e
//         sparsematrixvectormultiply( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                                     (const Float*)u, tempM, nShadowedU, shadowU );
//         scale( tempM, M, -1. );
//         add( tempM, (const Float*)e, M );
//         setfloats( M, sigma, 0. );
//         FREEZE( tempM, nShadowedSigma, shadowSigma );
//         ConjugateResidualSolverCSR( M, sigma, tempM, csrrowsA, csrcolumnsA, entriesA,
//                                     ressigma, nShadowedSigma, shadowSigma );
//         
//         // r = d = f - B sigma - C u
//         copyfloats( r, f, pdim );
//         sparsematrixvectormultiply( pdim, M, csrrowsB, csrcolumnsB, entriesB,
//                                     (const Float*)sigma, tempN, nShadowedSigma, shadowSigma );
//         addscaled( r, (const Float*)tempN, pdim, -1. );
//         sparsematrixvectormultiply( pdim, pdim, csrrowsC, csrcolumnsC, entriesC,
//                                     (const Float*)u, tempN, nShadowedU, shadowU );
//         addscaled( r, (const Float*)tempN, pdim, -1. );
//         FREEZE( r, nShadowedU, shadowU );
//         
// //         printf( "@@@@@@@@@@ Other residual sigma %f \n", vectornorm( M, ressigma ) );
// //         printf( "@@@@@@@@@@ Other residual u     %f \n", vectornorm( pdim, resu  ) );
//     
//     }
//     
//     if(false){
//     
//         copyfloats( ressigma, (const Float*)e, M );
//     
//         sparsematrixvectormultiply( M, M, csrrowsA, csrcolumnsA, entriesA,
//                                     (const Float*)sigma, tempM, nShadowedSigma, shadowSigma );
//                                     
//         sub( ressigma, (const Float*)tempM, M );
//         
//         sparsematrixvectormultiply( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                                     (const Float*)u, tempM, nShadowedU, shadowU );
//                                     
//         sub( ressigma, (const Float*)tempM, M );
//         
//         printf( "@@@@@@@@@@ One more Residual sigma %f \n", vectornorm( M, ressigma ) );
//         printf( "@@@@@@@@@@ One more Residual u     %f \n", vectornorm( pdim, resu  ) );
//     
//     }
//     
//     
//     
//     
//     
//     /******************/
//     /*** END **********/
//     /******************/
//     
//     printline( "@@@@@@@@@@ Algorithm finished" );
//     
//     free( r       );
//     free( C_r     ); 
//     free( AiBt_r  );
//     free( BAiBt_r );
//     free( d       );
//     free( C_d     );
//     free( AiBt_d  );
//     free( BAiBt_d );
//     free( p       );
//     free( C_p     );
//     free( AiBt_p  );
//     free( BAiBt_p );
//     free( tempM );
//     free( tempN );
//     free( entriesBt );
//     free( csrrowsBt );
//     free( csrcolumnsBt );
//     
//     printf( "@@@@@@@@@@ Exit: %f %f\n", vectornorm( M, ressigma ), vectornorm( pdim, resu  ) );
//     
//     
// }                            



// void UzawaConjugateResidualCSR_alt( const int M, const int pdim,
//                                 Float* const sigma, Float* const u,
//                                 Float* const ressigma, Float *const resu,
//                                 const Float* const e, const Float* const f,
//                                 int nShadowedSigma, int nShadowedU,
//                                 const int* shadowSigma, const int* shadowU,
//                                 const Float* const entriesA, const int* const csrrowsA, const int* const csrcolumnsA,
//                                 const Float* const entriesB, const int* const csrrowsB, const int* const csrcolumnsB,
//                                 const Float* const entriesC, const int* const csrrowsC, const int* const csrcolumnsC
//                             )
// {
//     
//     /**************************
//       A Bt | s  =  e
//       B C  | u  =  f
//     ***************************
//       exactly this sign
//       convention!!!
//     **************************/
//     assert( M > 0 );
//     assert( pdim > 0 );
//     assert( sigma );
//     assert( u );
//     assert( e );
//     assert( f );
//     assert( nShadowedSigma >= 0 );
//     assert( nShadowedU >= 0 );
//     if( nShadowedSigma > 0 ) assert( shadowSigma );
//     if( nShadowedU > 0 ) assert( shadowU );
//     assert( entriesA ); assert( csrrowsA ); assert( csrcolumnsA );
//     assert( entriesB ); assert( csrrowsB ); assert( csrcolumnsB );
//     assert( entriesC ); assert( csrrowsC ); assert( csrcolumnsC );
//     
//     Float* entriesBt = NULL;
//     int*   csrrowsBt = NULL;
//     int*   csrcolumnsBt = NULL;
//     
//     sparsetranspose( pdim, M,
//         (const int*)csrrowsB, (const int*)csrcolumnsB, (const Float*)entriesB,
//         &csrrowsBt, &csrcolumnsBt, &entriesBt );
//     
//     check_matrix( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt );
//     check_matrix( pdim, M, csrrowsB,  csrcolumnsB,  entriesB );
//     check_matrix( pdim, pdim, csrrowsC, csrcolumnsC, entriesC );
//     check_matrix( M, M, csrrowsA, csrcolumnsA, entriesA );
//     
//     Float* tempM = (Float*)mallocfloat(    M ); assert(    M );
//     Float* tempN = (Float*)mallocfloat( pdim ); assert( pdim );
//     
//     setfloats( M,    ressigma, 0. );
//     setfloats( pdim, resu,     0. );
//     setfloats( M,    tempM,    0. );
//     setfloats( pdim, tempN,    0. );
//     
//     Float* r       = (Float*)mallocfloat( pdim ); assert(       r );
//     Float* C_r     = (Float*)mallocfloat( pdim ); assert(     C_r );
//     Float* AiBt_r  = (Float*)mallocfloat(    M ); assert(  AiBt_r );
//     Float* BAiBt_r = (Float*)mallocfloat( pdim ); assert( BAiBt_r );
//     Float* d       = (Float*)mallocfloat( pdim ); assert(       d );
//     Float* C_d     = (Float*)mallocfloat( pdim ); assert(     C_d );
//     Float* AiBt_d  = (Float*)mallocfloat(    M ); assert(  AiBt_d );
//     Float* BAiBt_d = (Float*)mallocfloat( pdim ); assert( BAiBt_d );
//     Float* p       = (Float*)mallocfloat( pdim ); assert(       p );
//     Float* C_p     = (Float*)mallocfloat( pdim ); assert(     C_p );
//     Float* AiBt_p  = (Float*)mallocfloat(    M ); assert(  AiBt_p );
//     Float* BAiBt_p = (Float*)mallocfloat( pdim ); assert( BAiBt_p );
// 
//     setfloats( pdim,       r, 0. );
//     setfloats( pdim,     C_r, 0. );
//     setfloats(    M,  AiBt_r, 0. );
//     setfloats( pdim, BAiBt_r, 0. );
//     setfloats( pdim,       d, 0. );
//     setfloats( pdim,     C_d, 0. );
//     setfloats(    M,  AiBt_d, 0. );
//     setfloats( pdim, BAiBt_d, 0. );
//     setfloats( pdim,       p, 0. );
//     setfloats( pdim,     C_p, 0. );
//     setfloats(    M,  AiBt_p, 0. );
//     setfloats( pdim, BAiBt_p, 0. );
//     
//     int maxiter = 10 * (M + pdim);
//     int iter = 0;
//     
//     const Float allowed_error = allowed_error;
//     
//     /*************/
//     /* MAIN LOOP */
//     /*************/
//     
//     while( iter < maxiter ) {
//         
//         printf("@@@@@@@@@@ Uzawa-CRM Iteration %d / %d: unorm is %f\n", iter, maxiter, vectornorm(pdim,u) );
//         
//         /* Start or Restart condition check */
//         if( iter == 0 or ( iter % 1000 == 0 ) )
//         {
//             
//             // tempM = B^t u, tempM = -tempM, tempM = tempM + e
//             applyoperator( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                            (const Float*)u, tempM, nShadowedSigma, shadowSigma, nShadowedU, shadowU );
//             scale( tempM, M, -1. );
//             add( tempM, (const Float*)e, M );
//             // ...
//             
//             // sigma = A^\inv tempM
//             setfloats( M, sigma, 0. );
//             ConjugateResidualSolverCSR( M, sigma, tempM, csrrowsA, csrcolumnsA, entriesA,
//                                         ressigma, nShadowedSigma, shadowSigma );
//             
//             // r = d = f - B sigma - C u
//             copyfloats( r, f, pdim );
//             // ...
//             
//             applyoperator( pdim, M, csrrowsB, csrcolumnsB, entriesB,
//                            (const Float*)sigma, tempN, nShadowedU, shadowU, nShadowedSigma, shadowSigma );
//             addscaled( r, (const Float*)tempN, pdim, -1. );
//             applyoperator( pdim, pdim, csrrowsC, csrcolumnsC, entriesC,
//                            (const Float*)u, tempN, nShadowedU, shadowU, nShadowedU, shadowU );
//             addscaled( r, (const Float*)tempN, pdim, -1. );
//             
//             // residual Auxiliary variables
//             applyoperator( pdim, pdim, csrrowsC, csrcolumnsC, entriesC,
//                            (const Float*)r, C_r, nShadowedU, shadowU, nShadowedU, shadowU );
//             applyoperator( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                            (const Float*)r, tempM, nShadowedSigma, shadowSigma, nShadowedU, shadowU );
//             setfloats( M, AiBt_r, 0. );
//             ConjugateResidualSolverCSR( M, AiBt_r, tempM, csrrowsA, csrcolumnsA, entriesA,
//                                         ressigma, nShadowedSigma, shadowSigma );
//             applyoperator( pdim,    M,  csrrowsB,  csrcolumnsB,  entriesB,
//                            (const Float*)AiBt_r, BAiBt_r, nShadowedU, shadowU, nShadowedSigma, shadowSigma );
//             
//             // direction auxiliary variables
//             cpyfloats( pdim,       d, (const Float*)r       );
//             cpyfloats( pdim,     C_d, (const Float*)C_r     );
//             cpyfloats(    M,  AiBt_d, (const Float*)AiBt_r  );
//             cpyfloats( pdim, BAiBt_d, (const Float*)BAiBt_r );
//             
//         }
// 
//         Float residualnorm = vectornorm( pdim, r );
// 
//         printf( "@@@@@@@@@@ Residual norm: %f\n", residualnorm );
//         
//         if( residualnorm < allowed_error ) {
//             printf("@@@@@@@@@@ Threshold deceeded.\n");
//             break;
//         }
//         
//         
//         /* compute first variables */
//         
//         // t1
//         Float t1 = scalarproduct( pdim, r, C_r ) - scalarproduct( pdim, r, BAiBt_r );
//         
//         // p, C_p, AiBt_p, BAiBt_p == M d, C M d, AiBt M d, B Ai Bt M d
//         copyfloats( p, (const Float*)C_d,     pdim );
//         sub       ( p, (const Float*)BAiBt_d, pdim );
//         applyoperator( pdim, pdim,  csrrowsC,  csrcolumnsC,  entriesC,
//                        (const Float*)p, C_p, nShadowedU, shadowU, nShadowedU, shadowU );
//         applyoperator( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                        (const Float*)p, tempM, nShadowedSigma, shadowSigma, nShadowedU, shadowU );
//         
//         setfloats( M, AiBt_p, 0. );
//         ConjugateResidualSolverCSR( M, AiBt_p, (const Float*)tempM, csrrowsA, csrcolumnsA, entriesA,
//                                     ressigma, nShadowedSigma, shadowSigma );
//         applyoperator( pdim, M, csrrowsB, csrcolumnsB, entriesB,
//                        (const Float*)AiBt_p, BAiBt_p, nShadowedU, shadowU, nShadowedSigma, shadowSigma );
//         
//         // t2
//         Float t2 = scalarproduct( pdim, d, C_p ) - scalarproduct( pdim, d, BAiBt_p );
//         
//         // alpha
//         Float alpha = t1 / t2;
//         
//         // Update u, sigma
//         addscaled(       u,  (const Float*)     d, pdim,  alpha );
//         addscaled(   sigma,  (const Float*)AiBt_d,    M, -alpha );
//         
//         // Update residual terms
//         addscaled(       r, (const Float*)      p, pdim, -alpha );
//         addscaled(     C_r, (const Float*)    C_p, pdim, -alpha );
//         addscaled(  AiBt_r, (const Float*) AiBt_p,    M, -alpha );
//         addscaled( BAiBt_r, (const Float*)BAiBt_p, pdim, -alpha );
//         
//         
// /*        
//         
//         
//         if(false)
//         {
//             
//             // tempM = B^t u, tempM = -tempM, tempM = tempM + e
//             applyoperator( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                            (const Float*)u, tempM, nShadowedSigma, shadowSigma, nShadowedU, shadowU );
//             
//             scale( tempM, M, -1. );
//             
//             add( tempM, (const Float*)e, M );
//             // ...
//             
//             // sigma = A^\inv tempM
//             setfloats( M, sigma, 0. );
//             ConjugateResidualSolverCSR( M, sigma, tempM, csrrowsA, csrcolumnsA, entriesA,
//                                         ressigma, nShadowedSigma, shadowSigma );
//             
//             // r = d = f - B sigma - C u
//             copyfloats( r, f, pdim );
//             // ...
//             
//             applyoperator( pdim, M, csrrowsB, csrcolumnsB, entriesB,
//                            (const Float*)sigma, tempN, nShadowedU, shadowU, nShadowedSigma, shadowSigma );
//                                         
//             addscaled( r, (const Float*)tempN, pdim, -1. );
//             
//             applyoperator( pdim, pdim, csrrowsC, csrcolumnsC, entriesC,
//                            (const Float*)u, tempN, nShadowedU, shadowU, nShadowedU, shadowU );
//                                         
//             addscaled( r, (const Float*)tempN, pdim, -1. );
//             
//         }*/
//         
// //         if(false){
// //             
// //             // residual Auxiliary variables
// //             applyoperator( pdim, pdim, csrrowsC, csrcolumnsC, entriesC,
// //                            (const Float*)r, C_r, nShadowedU, shadowU, nShadowedU, shadowU );
// //                                              
// //             applyoperator( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
// //                            (const Float*)r, tempM, nShadowedSigma, shadowSigma, nShadowedU, shadowU );
// //             
// //             setfloats( M, AiBt_r, 0. );
// //             ConjugateResidualSolverCSR( M, AiBt_r, tempM, csrrowsA, csrcolumnsA, entriesA,
// //                                         ressigma, nShadowedSigma, shadowSigma );
// //                                         
// //             applyoperator( pdim,    M,  csrrowsB,  csrcolumnsB,  entriesB,
// //                            (const Float*)AiBt_r, BAiBt_r, nShadowedU, shadowU, nShadowedSigma, shadowSigma );
//             
// //             // direction auxiliary variables
// //             cpyfloats( pdim,       d, (const Float*)r       );
// //             cpyfloats( pdim,     C_d, (const Float*)C_r     );
// //             cpyfloats(    M,  AiBt_d, (const Float*)AiBt_r  );
// //             cpyfloats( pdim, BAiBt_d, (const Float*)BAiBt_r );
// //            
// //        }
//         
//         
//         
//         
//         
//         
//         // Compute beta
//         Float beta = scalarproduct( pdim, r, C_r ) - scalarproduct( pdim, r, BAiBt_r ); 
//         
//         // Update direction terms
//         scalefloats( pdim,       d, beta );
//         scalefloats( pdim,     C_d, beta );
//         scalefloats(    M,  AiBt_d, beta );
//         scalefloats( pdim, BAiBt_d, beta );
//         
//         add(       d, (const Float*)      r, pdim );
//         add(     C_d, (const Float*)    C_r, pdim );
//         add(  AiBt_d, (const Float*) AiBt_r,    M );
//         add( BAiBt_d, (const Float*)BAiBt_r, pdim );
//         
//         iter++;
//         
//         // direction auxiliary variables
//         cpyfloats( pdim,       d, (const Float*)r       );
//         cpyfloats( pdim,     C_d, (const Float*)C_r     );
//         cpyfloats(    M,  AiBt_d, (const Float*)AiBt_r  );
//         cpyfloats( pdim, BAiBt_d, (const Float*)BAiBt_r );
//         
//     }
//     
//     /******************/
//     /*** END LOOP *****/
//     /******************/
//     
//     {
//             
//         copyfloats( ressigma, (const Float*)e, M );
//         FREEZE( ressigma, nShadowedSigma, shadowSigma );
//     
//         FREEZE( tempM, nShadowedSigma, shadowSigma );
//         applyoperator( M, M, csrrowsA, csrcolumnsA, entriesA,
//                        (const Float*)sigma, tempM, nShadowedSigma, shadowSigma, nShadowedSigma, shadowSigma );
//                                     
//         FREEZE( sigma, nShadowedSigma, shadowSigma );
//         sub( ressigma, (const Float*)tempM, M );
//         
//         FREEZE( u, nShadowedU, shadowU );
//         applyoperator( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                        (const Float*)u, tempM, nShadowedSigma, shadowSigma, nShadowedU, shadowU );
//                                     
//         sub( ressigma, (const Float*)tempM, M );
//         FREEZE( ressigma, nShadowedSigma, shadowSigma );
//         
//         copyfloats( resu, r, pdim );
//         FREEZE( resu, nShadowedU, shadowU );
//         
//         
//     }
//     
// //     printf( "@@@@@@@@@@ Residual sigma %f \n", vectornorm( M, ressigma ) );
// //     printf( "@@@@@@@@@@ Residual u     %f \n", vectornorm( pdim, resu  ) );
//     
//     if(true){
//     
//         // sigma and ressigma
//         // tempM = B^t u, tempM = -tempM, tempM = tempM + e
//         applyoperator( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                        (const Float*)u, tempM, nShadowedSigma, shadowSigma, nShadowedU, shadowU );
//         scale( tempM, M, -1. );
//         add( tempM, (const Float*)e, M );
//         setfloats( M, sigma, 0. );
//         FREEZE( tempM, nShadowedSigma, shadowSigma );
//         ConjugateGradientSolverCSR_Formatted( M, sigma, tempM, csrrowsA, csrcolumnsA, entriesA,
//                                     ressigma, nShadowedSigma, shadowSigma );
//         
//         // r = d = f - B sigma - C u
//         copyfloats( r, f, pdim );
//         applyoperator( pdim, M, csrrowsB, csrcolumnsB, entriesB,
//                        (const Float*)sigma, tempN, nShadowedU, shadowU, nShadowedSigma, shadowSigma );
//         addscaled( r, (const Float*)tempN, pdim, -1. );
//         applyoperator( pdim, pdim, csrrowsC, csrcolumnsC, entriesC,
//                        (const Float*)u, tempN, nShadowedU, shadowU, nShadowedU, shadowU );
//         addscaled( r, (const Float*)tempN, pdim, -1. );
//         FREEZE( r, nShadowedU, shadowU );
//         
// //         printf( "@@@@@@@@@@ Other residual sigma %f \n", vectornorm( M, ressigma ) );
// //         printf( "@@@@@@@@@@ Other residual u     %f \n", vectornorm( pdim, resu  ) );
//     
//     }
//     
//     if(false){
//     
//         copyfloats( ressigma, (const Float*)e, M );
//     
//         applyoperator( M, M, csrrowsA, csrcolumnsA, entriesA,
//                        (const Float*)sigma, tempM, nShadowedSigma, shadowSigma, nShadowedSigma, shadowSigma );
//                                     
//         sub( ressigma, (const Float*)tempM, M );
//         
//         applyoperator( M, pdim, csrrowsBt, csrcolumnsBt, entriesBt,
//                        (const Float*)u, tempM, nShadowedSigma, shadowSigma, nShadowedU, shadowU );
//                                     
//         sub( ressigma, (const Float*)tempM, M );
//         
//         printf( "@@@@@@@@@@ One more Residual sigma %f \n", vectornorm( M, ressigma ) );
//         printf( "@@@@@@@@@@ One more Residual u     %f \n", vectornorm( pdim, resu  ) );
//     
//     }
//     
//     
//     
//     
//     
//     /******************/
//     /*** END **********/
//     /******************/
//     
//     printline( "@@@@@@@@@@ Algorithm finished" );
//     
//     free( r       );
//     free( C_r     ); 
//     free( AiBt_r  );
//     free( BAiBt_r );
//     free( d       );
//     free( C_d     );
//     free( AiBt_d  );
//     free( BAiBt_d );
//     free( p       );
//     free( C_p     );
//     free( AiBt_p  );
//     free( BAiBt_p );
//     free( tempM );
//     free( tempN );
//     free( entriesBt );
//     free( csrrowsBt );
//     free( csrcolumnsBt );
//     
//     printf( "@@@@@@@@@@ Exit: %f %f\n", vectornorm( M, ressigma ), vectornorm( pdim, resu  ) );
//     
//     
// }         



