#ifndef INCLUDEGUARD_SOLVER_SYSTEMSOLVER
#define INCLUDEGUARD_SOLVER_SYSTEMSOLVER

#include "../basic.hpp"
#include "sparsesolver.hpp"

void HodgeConjugateResidualSolverCSR( 
    const int N, 
    const int L, 
    Float* x, 
    const Float* b, 
    const int*  Arows, const int*  Acolumns, const Float*  Avalues, 
    const int*  Brows, const int*  Bcolumns, const Float*  Bvalues, 
    const int* Btrows, const int* Btcolumns, const Float* Btvalues, 
    Float* residual,
    Float allowed_error,
    unsigned int restart_modulo
);






#include <cmath>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <omp.h>
#include <utility>




 

void HodgeConjugateResidualSolverCSR( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    Float* res,
    Float allowed_error,
    unsigned int restart_modulo
) {
    
    assert( N > 0 );
    assert( L > 0 );
    assert( x );
    assert( b );
    assert(  Arows );
    assert(  Acolumns );
    assert(  Avalues );
    assert(  Brows );
    assert(  Bcolumns );
    assert(  Bvalues );
    assert( Btrows );
    assert( Btcolumns );
    assert( Btvalues );
    assert( res );
    assert( allowed_error > 0 );
    assert( restart_modulo > 0 );
    
    Float* __restrict__  dir = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ Mdir = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ Mres = (Float*)malloc( sizeof(Float) * N );
    
    Float* __restrict__ aux1 = (Float*)malloc( sizeof(Float) * L );
    Float* __restrict__ aux2 = (Float*)malloc( sizeof(Float) * L );
    Float* __restrict__ auxR = (Float*)malloc( sizeof(Float) * L );
    
    Float* __restrict__  vil = (Float*)malloc( sizeof(Float) * N );
    
    assert(  dir );
    assert( Mdir );
    assert( Mres );
    
    assert( aux1 );
    assert( aux2 );
    assert( auxR );
    
    assert( vil );
    
    
    Float Md_r;
    Float Md_Md;
    
    int k = 0;
    
    while( k < N ){
        
        bool restart_condition = ( k == 0 ); // or k % 1000 == 0;
        
        bool r_seems_small = false;//std::sqrt(Md_r) < allowed_error;

        if( restart_condition or r_seems_small ) {
            
            #pragma omp parallel for
            for( int c = 0; c < L; c++ ) {
                
                aux1[c] = 0.;
                
                for( int d = Btrows[c]; d < Btrows[c+1]; d++ )
                    aux1[c] += Btvalues[ d ] * x[ Btcolumns[d] ];
                
            }
            
            for( int c = 0; c < L; c++ ) aux2[c] = 0.;
            ConjugateGradientSolverCSR( 
                L, 
                aux2, 
                (const Float *)aux1, 
                Arows, Acolumns, Avalues, 
                auxR,
                allowed_error,
                restart_modulo
            );
            
            #pragma omp parallel for
            for( int c = 0; c < N; c++ ) {
                
                res[c] = b[c];
                
                for( int d = Brows[c]; d < Brows[c+1]; d++ )
                    res[c] -= Bvalues[ d ] * aux2[ Bcolumns[d] ];
                
                dir[c] = res[c]; 
            
            }
            
            Md_r  = 0.;
            Md_Md = 0.;
            
            #pragma omp parallel for
            for( int c = 0; c < L; c++ ) {
                
                aux1[c] = 0;
                
                for( int d = Btrows[c]; d < Btrows[c+1]; d++ )
                    aux1[c] += Btvalues[ d ] * dir[ Btcolumns[d] ];
                
            }
            
            for( int c = 0; c < L; c++ ) aux2[c] = 0.;
            ConjugateGradientSolverCSR( 
                L, 
                aux2, 
                (const Float *)aux1, 
                Arows, Acolumns, Avalues, 
                auxR,
                allowed_error,
                restart_modulo
            );
            
            #pragma omp parallel for reduction(+:Md_r,Md_Md)
            for( int c = 0; c < N; c++ ) {
                
                Mdir[c] = 0.;
                
                for( int d = Brows[c]; d < Brows[c+1]; d++ )
                    Mdir[c] += Bvalues[ d ] * aux2[ Bcolumns[d] ];
                
                Mres[c] = Mdir[c];
                
                Md_r  += Mdir[c] *  res[c];
                Md_Md += Mdir[c] * Mdir[c];
                
            }
            
            
        }
        
        /* Check whether res is small */
                
        bool r_is_small = std::sqrt(Md_r) < allowed_error;
        
        if( r_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        Float alpha = Md_r / Md_Md;
        
        Float new_Mr_r = 0.;
        
        #pragma omp parallel for
        for( int c = 0; c < L; c++ ) {
            
            aux1[c] = 0;
            
            for( int d = Btrows[c]; d < Btrows[c+1]; d++ )
                aux1[c] += Btvalues[ d ] * Mdir[ Btcolumns[d] ];
            
        }
        
        // for( int c = 0; c < L; c++ ) aux2[c] = 0.;
        ConjugateGradientSolverCSR( 
            L, 
            aux2, 
            (const Float *)aux1, 
            Arows, Acolumns, Avalues, 
            auxR,
            allowed_error,
            restart_modulo
        );
        
        #pragma omp parallel for reduction(+:new_Mr_r)
        // #pragma omp parallel for 
        for( int c = 0; c < N; c++ ) {
            
            vil[c] = 0.;
            
            for( int d = Brows[c]; d < Brows[c+1]; d++ )
                vil[c] += Bvalues[ d ] * aux2[ Bcolumns[d] ];
            
        // }
                    
        // #pragma omp parallel for reduction(+:new_Ar_r)
        // for( int c = 0; c < N; c++ )
        // {
        
            x[c]    =   x[c] + alpha *  dir[c];
            
            res[c]  = res[c] - alpha * Mdir[c];
            
            Mres[c] = Mres[c] - alpha * vil[c];
            
            new_Mr_r += Mres[c] * res[c];
            
        }
        
        Float beta = new_Mr_r / Md_r;
        
        Md_Md = 0.;
        #pragma omp parallel for reduction(+:Md_Md)
        for( int c = 0; c < N; c++ )
        {
            
             dir[c] =  res[c] + beta *  dir[c];
            
            Mdir[c] = Mres[c] + beta * Mdir[c];
            
            Md_Md += Mdir[c] * Mdir[c];
        }
        
        Md_r = new_Mr_r;
        
        printf("Hodge Residual after %d of max. %d iterations: %.9Le\n", k, N, (long double)std::sqrt(Md_r) );
        
//         if( k % 100 == 0 ) 
//         printf("At Iteration %d we have %.9Le --- [%.9Le,%.9Le,%.9Le,%.9Le]\n",
//             k,
//             (long double)std::sqrt(r_r), (long double)alpha, (long double)beta,
//             (long double)d_Ad, (long double)r_r_new );
        
        k++;
        
    };
    
    printf("Hodge Residual after %d of max. %d iterations: %.9Le\n", k, N, (long double)std::sqrt(Md_r) );

    
    free(  dir );
    free( Mdir );
    free( Mres );

    free( aux1 );
    free( aux2 );
    free( auxR );
    
    free( vil );

}
 
 
 

 
#endif
