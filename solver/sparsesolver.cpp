


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
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float allowed_error
) {
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    Float* direction = (Float*)malloc( sizeof(Float) * N );
    Float* auxiliary = (Float*)malloc( sizeof(Float) * N );
    assert( direction );
    assert( auxiliary );
    

    #pragma omp parallel for
    for( int c = 0; c < N; c++ ) {
        direction[c] = 0.;
        residual[c] = 0;
        auxiliary[c] = 0.;
    }
    
    #pragma omp parallel for
    for( int c = 0; c < N; c++ ) {

        for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
            residual[c] = b[c] - 1. * csrvalues[ d ] * x[ csrcolumns[d] ];
        
        direction[c] = residual[c];
    }
    

    const Float eps = allowed_error;
    int K = 0;
    int emergency = 0;
    
    while( K < N ){
        K++;
        
        Float gamma = 0.;
        Float alpha = 0.;
        Float beta  = 0.;
        Float temp1 = 0.;
        Float temp2 = 0.;
        
        /* Compute residual norm */
        /* After some iterations, restart descent*/
        if( emergency || (K % 1000 == 0) ) {
    
            emergency = 0;
    
            #pragma omp parallel for reduction(+:gamma)
            for( int c = 0; c < N; c++ ) {
                
                residual[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                gamma += residual[c] * residual[c];
            }
            
        } else {
            
            #pragma omp parallel for reduction(+:gamma)
            for( int c = 0; c < N; c++ )
                gamma += residual[c] * residual[c];
                
        }
        
        gamma = sqrt(gamma);
        
        
        /* Check whether residual is small */
        if( gamma < eps ) {

            for( int c = 0; c < N; c++ ) {
                residual[c] = b[c];
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
            }

            Float residualnorm = 0.;

            for( int c = 0; c < N; c++ )
                residualnorm += residual[c]*residual[c];
            
            break;

        };
        
        #pragma omp parallel for reduction(+:alpha,temp1)
        for( int c = 0; c < N; c++ )
        {
            alpha += direction[c] * residual[c];
            
            auxiliary[c] = 0.;
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                auxiliary[c] += csrvalues[ d ] * direction[ csrcolumns[d] ];
                    
            temp1 += direction[c] * auxiliary[c];
        }
        
        alpha /= temp1;
        
        if( temp1 < 0.00000001 ) emergency++;
                
        #pragma omp parallel for reduction(+:temp2,beta)
        for( int c = 0; c < N; c++ )
        {
            
            x[c] += alpha * direction[c];
            
            temp2 += residual[c] * residual[c];
            
            residual[c] -= alpha * auxiliary[c];
            
            beta += residual[c] * residual[c];
        }
        
        beta /= temp2;
        
        if( temp2 < 0.00000001 ) emergency++;
        

        #pragma omp parallel for
        for( int c = 0; c < N; c++ )
            auxiliary[c] = residual[c] + beta * direction[c];
        
        
                
        {
            Float* swappi = auxiliary;
            auxiliary = direction;
            direction = swappi;
        }
        
        printf("At Iteration %d we have %Lf --- [%Lf,%Lf,%Lf,%Lf]\n", K,
                     (long double)gamma, (long double)alpha, (long double)beta,
                     (long double)temp1, (long double)temp2 );
        
        if( emergency ) {
            puts("Emergency!");
        }
        
    };
    
    for( int c = 0; c < N; c++ ) {
        residual[c] = b[c];
        for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
            residual[c] += -1. * csrvalues[ d ] * x[ csrcolumns[d] ];
    }

    Float residualnorm = 0.;
    for( int c = 0; c < N; c++ ) {
        
        residual[c] = b[c];
        for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
            residual[c] += -1. * csrvalues[ d ] * x[ csrcolumns[d] ];
        residualnorm += residual[c]*residual[c];
    }
    printf("Residual after %d of max. %d iterations: %Lf\n", K, N, (long double)sqrt(residualnorm) );

    
    free( direction ); 
    free( auxiliary );

}




void ConjugateResidualSolverCSR( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float allowed_error
) {
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    Float* direction = (Float*)malloc( sizeof(Float) * N );
    Float* auxiliary = (Float*)malloc( sizeof(Float) * N );
    assert( direction );
    assert( auxiliary );
    

    #pragma omp parallel for
    for( int c = 0; c < N; c++ ) {
        direction[c] = 0.;
        residual[c] = 0.;
        auxiliary[c] = 0.;
    }
    
    #pragma omp parallel for
    for( int c = 0; c < N; c++ ) {
        
        residual[c] = b[c];
        for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
        {
            assert( 0 <= csrcolumns[d] );
            assert( csrcolumns[d] <= N );
                residual[c] += -1. * csrvalues[ d ] * x[ csrcolumns[d] ];
        }
        direction[c] = residual[c];
    }
    
    {

        Float residualnorm = 0.;
        for( int c = 0; c < N; c++ ) {
            residual[c] = b[c];
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                residual[c] += -1. * csrvalues[ d ] * x[ csrcolumns[d] ];
            residualnorm += residual[c]*residual[c];
        }
    
    }
    
    const Float eps = allowed_error;
    int K = 0;
    int emergency = 0;
    
    while( K < N ){
        K++;
        
        Float gamma = 0.;
        
        /* Compute residual norm */
        /* After some iterations, restart descent*/
        if( emergency || (K % 1000 == 0) || (K > N - 10 ) ) {
            
            emergency = 0;
            
            #pragma omp parallel for reduction(+:gamma)
            for( int c = 0; c < N; c++ ) {
                residual[c] = b[c];
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                gamma += residual[c]*residual[c];
            }
            
        } else {
        
            #pragma omp parallel for reduction(+:gamma)
            for( int c = 0; c < N; c++ )
                gamma += residual[c] * residual[c];
        
        }
        
        gamma = sqrt(gamma);
        
        
        
        /* Check whether residual is small */
        if( gamma < eps ) {

            for( int c = 0; c < N; c++ ) {
                
                residual[c] = b[c];
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
            }
            break;
        };
        
        
        Float residualanorm = 0.;
        for( int c = 0; c < N; c++ ) {
            
            Float help = 0.;
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                help += csrvalues[ d ] * direction[ csrcolumns[d] ];
            residualanorm += residual[c] * help;
        }
        
        Float temp1 = 0.;
        for( int c = 0; c < N; c++ ) {
            
            Float help = 0.;
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                help += csrvalues[ d ] * direction[ csrcolumns[d] ];
            temp1 += help * help;
        }
        
        Float alpha = residualanorm / temp1;
        
        for( int c = 0; c < N; c++ ) {
            
            x[c] += alpha * direction[c];
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                residual[c] -= alpha * csrvalues[ d ] * direction[ csrcolumns[d] ];
        }
        
        Float residualanorm2 = 0.;
        for( int c = 0; c < N; c++ ) {
            
            Float help = 0.;
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                help += csrvalues[ d ] * residual[ csrcolumns[d] ];
            residualanorm2 += residual[c] * help;
        }
        
        Float beta = residualanorm2 / residualanorm;
        
        for( int c = 0; c < N; c++ ) {
            
            direction[c] = residual[c] + beta * direction[c];
        }
        
        
        if( true )
            printf("At Iteration %d we have %Lf --- [%Lf,%Lf,%Lf,%Lf]\n", K,
                     (long double)gamma, (long double)alpha, (long double)beta,
                     (long double)residualanorm, (long double)residualanorm2 );
        
        if( emergency ) {
//             puts("Emergency!");
        }
        
        if( fabs( residualanorm2 ) < eps ) break;
        
    };
    
    Float residualnorm = 0.;
    for( int c = 0; c < N; c++ ) {
        
        residual[c] = b[c];
        for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
            residual[c] += -1. * csrvalues[ d ] * x[ csrcolumns[d] ];
        residualnorm += residual[c]*residual[c];
    }
    printf("Residual after %d of max. %d iterations: %Lf\n", K, N, (long double)sqrt(residualnorm) );
    

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
//     const Float eps = allowed_error;
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
//         if( iter == 0 || ( false && iter % 1000 == 0 ) ) {
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
//         if( residualnorm < eps ) {
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
//     const Float eps = allowed_error;
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
//         if( iter == 0 || ( iter % 1000 == 0 ) )
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
//         if( residualnorm < eps ) {
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



