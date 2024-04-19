

// This one is taken from Wikipedia, I don't know what it is.

int WHATEVER( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ res,
    const Float tolerance,
    int print_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( res );
    assert( tolerance > 0 );
    assert( print_modulo >= -1 );
    
    Float* __restrict__  r = new (std::nothrow) Float[N];
    Float* __restrict__ p0 = new (std::nothrow) Float[N];
    Float* __restrict__ p1 = new (std::nothrow) Float[N];
    Float* __restrict__ p2 = new (std::nothrow) Float[N];
    Float* __restrict__ s0 = new (std::nothrow) Float[N];
    Float* __restrict__ s1 = new (std::nothrow) Float[N];
    Float* __restrict__ s2 = new (std::nothrow) Float[N];
    
    assert(  r );
    assert( p0 );
    assert( p1 );
    assert( p2 );
    assert( s0 );
    assert( s1 );
    assert( s2 );

    if( print_modulo >= 0 ) 
        LOGPRINTF( "START Whatever\n" );

    Float r_r = 0.;
    
    #if defined(_OPENMP)
    #pragma omp parallel for reduction(+:r_r)
    #endif 
    for( int c = 0; c < N; c++ ) {
        
        r[c] = b[c];
        
        for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
            r[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
        
        p0[c] =  r[c];
        p1[c] = p0[c];
        
        r_r += r[c] * r[c];
                
    }
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif 
    for( int c = 0; c < N; c++ ) {
        
        s0[c] = 0.;
        
        for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
            s0[c] += csrvalues[ d ] * p0[ csrcolumns[d] ];
        
        s1[c] = s0[c];
                
    }
    
    int K = 0;

    while( K < N ){
        
        Float  r_s1 = 0.;
        Float s0_s0 = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:r_s1,s0_s0)
        #endif 
        for( int c = 0; c < N; c++ ) {
            
            p2[c] = p1[c];
            p1[c] = p0[c];
            
            s2[c] = s1[c];
            s1[c] = s0[c];
            
             r_s1 +=  r[c] * s1[c];
            s0_s0 += s1[c] * s1[c];
                    
        }
        
        Float alpha = r_s1 / s0_s0;
    
        r_r = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:r_r)
        #endif 
        for( int c = 0; c < N; c++ ) {
            
            x[c] += alpha * p1[c];
            r[c] -= alpha * s1[c];
            
            r_r += r[c] * r[c];
        }
        
        if( std::sqrt(r_r) < tolerance )
            break;
        
        Float s0_s1 = 0.;
        Float s1_s1 = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:s0_s1,s1_s1)
        #endif 
        for( int c = 0; c < N; c++ ) {
            
            p0[c] = s1[c];
            
            s0[c] = 0.;
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                s0[c] += csrvalues[ d ] * s1[ csrcolumns[d] ];
            
            s0_s1 += s0[c] * s1[c];
            s1_s1 += s1[c] * s1[c];
            
        }
        
        Float beta = s0_s1 / s1_s1;
        
        Float s0_s2 = 0.;
        Float s2_s2 = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:s0_s2,s2_s2)
        #endif 
        for( int c = 0; c < N; c++ ) {
            
            p0[c] -= beta * p1[c];
            s0[c] -= beta * s1[c];
            
            s0_s2 += s0[c] * s2[c];
            s2_s2 += s2[c] * s2[c];
        }
        
        if( K > 0 )
        {
            
            Float gamma = s0_s2 / s2_s2;
            
            #if defined(_OPENMP)
            #pragma omp parallel for 
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                p0[c] -= gamma * p2[c];
                s0[c] -= gamma * s2[c];
                
            }
            
        }
            
            
        
        if( print_modulo > 0 and K % print_modulo == 0 ) UNLIKELY 
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9Le < %.9Le\n", K, N, (long double)std::sqrt(r_r), (long double) tolerance );
        
        K++;
        
    }
    
    if( print_modulo >= 0 ) 
        LOGPRINTF( "(%d/%d)  FINISHED: Residual norm is %.9Le < %.9Le\n", K, N, (long double)std::sqrt(r_r), (long double) tolerance );

    
    delete[] (  r );
    delete[] ( p0 );
    delete[] ( p1 );
    delete[] ( p2 );
    delete[] ( s0 );
    delete[] ( s1 );
    delete[] ( s2 );

    return K;
    
}


