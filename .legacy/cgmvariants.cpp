

void ConjugateGradientSolverCSR_variant( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float allowed_error,
    unsigned int restart_modulo
);



void ConjugateGradientSolverCSR_sstep( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float allowed_error,
    unsigned int restart_modulo
);





















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
    
    Float*  h = (Float*)malloc( sizeof(Float) * N );
    Float*  d = (Float*)malloc( sizeof(Float) * N );
    Float* Ad = (Float*)malloc( sizeof(Float) * N );
    assert(  h );
    assert(  d );
    assert( Ad );
    
    Float alpha;
    Float  beta;
    Float   r_r;

    
    // Initialize : compute the residual, the direction, and r * r
    
    int K = 0;

    while( K < N ){
        
        
        bool restart_condition = (K == 0); // or (K % 1000 == 1);
        
        bool residual_seems_small = false and (std::sqrt(r_r) < allowed_error);
        
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
        
        
        #pragma omp parallel for 
        for( int c = 0; c < N; c++ )
        {
            
            h[c] = r[c] + beta * d[c];
            
        }
            

        Float Ad_d = 0.;
        
        #pragma omp parallel for reduction(+:Ad_d)
        for( int c = 0; c < N; c++ )
        {
            
//             h[c] = r[c] + beta * d[c];
            
            Ad[c] = 0.;
        
            for( int e = csrrows[c]; e < csrrows[c+1]; e++ )
                Ad[c] += csrvalues[ e ] * ( h[ csrcolumns[e] ] ); //= r[ csrcolumns[e] ] + beta * d[ csrcolumns[e] ] );
            
            Ad_d += h[c] * Ad[c];
            
            
        }
        
        alpha = r_r / Ad_d;
            
        Float r_r_new = 0.;
        
        #pragma omp parallel for reduction(+:r_r_new)
        for( int c = 0; c < N; c++ ) {
        
            x[c] = x[c] + alpha *  h[c];
            
            r[c] = r[c] - alpha * Ad[c];
            
            r_r_new += r[c] * r[c];
            
        }
        
        beta = r_r_new / r_r;
        
        r_r = r_r_new;
        
        std::swap(d,h);
        
//         if( K % 100 == 0 ) 
//         printf("At Iteration %d we have %.9Le --- [%.9Le,%.9Le,%.9Le,%.9Le]\n",
//             K,
//             (long double)std::sqrt(r_r), (long double)alpha, (long double)beta,
//             (long double)Ad_d, (long double)Ad_Ad );
        
        K++;
        
        
    };
    
    printf("Residual after %d of max. %d iterations: %.9Le\n", K, N, (long double)std::sqrt(r_r) );

    
    free(  h ); 
    free(  d ); 
    free( Ad );

}












// s-step iterative methods for symmetric linear systems. Chronopoulos & Gear



void ConjugateGradientSolverCSR_variant2( 
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

void ConjugateGradientSolverCSR_sstep_oneloop( 
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
    
    Float* __restrict__ d = (Float*)malloc( sizeof(Float) * N );
    Float*  r = (Float*)malloc( sizeof(Float) * N );
    Float*  Ad = (Float*)malloc( sizeof(Float) * N );
    Float*  Ar = (Float*)malloc( sizeof(Float) * N );
    assert(  d );
    assert(  r );
    assert( Ad );
    assert( Ar );
    
    Float*   r_new = (Float*)malloc( sizeof(Float) * N );
    Float*  Ad_new = (Float*)malloc( sizeof(Float) * N );
    Float*  Ar_new = (Float*)malloc( sizeof(Float) * N );
    assert(  r_new );
    assert( Ad_new );
    assert( Ar_new );
    
    
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
        
        Float Ar_r_new = 0.;
        
        #pragma omp parallel for reduction(+:r_r_new,Ar_r_new)
        for( int c = 0; c < N; c++ )
        {
            
            // compute the step
            
            d[c] = r[c] + beta * d[c];
            
            Ad_new[c] = Ar[c] + beta * Ad[c];
            
            // go along that step 
            
            x[c] = x[c] + alpha * d[c];
            
            r_new[c] = r[c] - alpha * Ad_new[c];
            
            r_r_new += r_new[c] * r_new[c];
            
            Ar_new[c] = 0.;
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                Ar_new[c] += csrvalues[ d ] * ( r[ csrcolumns[d] ] - alpha * ( Ar[ csrcolumns[d] ] + beta * Ad[ csrcolumns[d] ] ) );
            
            Ar_r_new += r_new[c] * Ar_new[c];
            
        }
        
        std::swap(  r,  r_new );
        std::swap( Ar, Ar_new );
        std::swap( Ad, Ad_new );
        
        beta = r_r_new / r_r;
        
        alpha = r_r_new / ( Ar_r_new - r_r_new * beta / alpha );
        
        r_r = r_r_new;
                
        
//         if( K % 100 == 0 ) 
//         printf("At Iteration %d we have %.9Le --- [%.9Le,%.9Le,%.9Le,%.9Le]\n",
//             K,
//             (long double)std::sqrt(r_r), (long double)alpha, (long double)beta,
//             (long double)Ad_d, (long double)d_r );
        
        K++;
        
        
    };
    
    printf("Residual after %d of max. %d iterations: %.9Le\n", K, N, (long double)std::sqrt(r_r) );

    for( int c = 0; c < N; c++ ) 
    {
        residual[c] =  r[c];
    }
    
    free(  d ); 
    free(  r ); 
    free( Ad );
    free( Ar );

    free(  r_new ); 
    free( Ad_new );
    free( Ar_new );

}





























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





















