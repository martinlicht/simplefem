

void ConjugateResidualSolverCSR_variant( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float allowed_error,
    unsigned int restart_modulo
);

void ConjugateResidualSolverCSR_old( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float allowed_error,
    unsigned int restart_modulo
);






void ConjugateResidualSolverCSR_variant( 
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























