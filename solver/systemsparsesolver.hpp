#ifndef INCLUDEGUARD_SOLVER_SYSTEMSOLVER
#define INCLUDEGUARD_SOLVER_SYSTEMSOLVER

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <utility>

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
    Float threshold,
    int print_modulo,
    Float inneriteration_threshold,
    int inneriteration_print_modulo
);

void HodgeConjugateResidualSolverCSR_textbook( 
    const int N, 
    const int L, 
    Float* x, 
    const Float* b, 
    const int*  Arows, const int*  Acolumns, const Float*  Avalues, 
    const int*  Brows, const int*  Bcolumns, const Float*  Bvalues, 
    const int* Btrows, const int* Btcolumns, const Float* Btvalues, 
    Float* residual,
    Float threshold,
    int print_modulo,
    Float inneriteration_threshold,
    int inneriteration_print_modulo
);

void HodgeConjugateResidualSolverCSR_SSOR( 
    const int N, 
    const int L, 
    Float* x, 
    const Float* b, 
    const int*  Arows, const int*  Acolumns, const Float*  Avalues, 
    const int*  Brows, const int*  Bcolumns, const Float*  Bvalues, 
    const int* Btrows, const int* Btcolumns, const Float* Btvalues, 
    Float* residual,
    Float threshold,
    int print_modulo,
    Float inneriteration_threshold,
    int inneriteration_print_modulo
);











 

void HodgeConjugateResidualSolverCSR( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    Float* res,
    Float threshold,
    int print_modulo,
    Float inneriteration_threshold,
    int inneriteration_print_modulo
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
    assert( threshold > 0 );
//     assert( print_modulo >= 0 );
    
    Float* __restrict__  dir = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ Mdir = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ Mres = (Float*)malloc( sizeof(Float) * N );
    
    Float* __restrict__ aux1 = (Float*)malloc( sizeof(Float) * L );
    Float* __restrict__ aux2 = (Float*)malloc( sizeof(Float) * L );
    Float* __restrict__ auxR = (Float*)malloc( sizeof(Float) * L );
    
    Float* __restrict__  vil = (Float*)malloc( sizeof(Float) * N );
    
    Float* __restrict__  precon = (Float*)malloc( sizeof(Float) * L );
    
    assert(  dir );
    assert( Mdir );
    assert( Mres );
    
    assert( aux1 );
    assert( aux2 );
    assert( auxR );
    
    assert( vil );
    
    assert( precon );
    
    
    #pragma omp parallel for
    for( int c = 0; c < L; c++ ) {
        
        precon[c] = 0.;
        
        for( int d = Arows[c]; d < Arows[c+1]; d++ )
            if( Acolumns[d] == c ) 
                precon[c] += Avalues[ d ];
            
        precon[c] = 1./precon[c];
//         precon[c] = 1.;
        
        assert( precon[c] > 0. );
        
    }
    
    
    
    Float Md_r  = notanumber;
    Float Md_Md = notanumber;
    
    int k = 0;
    
    while( k < N ){
        
        bool restart_condition = ( k == 0 ); // or k % 1000 == 0;
        
        bool residualenergy_seems_small = false;//absolute(Md_r) < threshold*threshold;

        if( restart_condition or residualenergy_seems_small ) {
            
            #pragma omp parallel for
            for( int c = 0; c < L; c++ ) {
                
                aux1[c] = 0.;
                
                for( int d = Btrows[c]; d < Btrows[c+1]; d++ )
                    aux1[c] += Btvalues[ d ] * x[ Btcolumns[d] ];
                
            }
            
            for( int c = 0; c < L; c++ ) aux2[c] = 0.;
            ConjugateGradientSolverCSR_DiagonalPreconditioner( 
                L, 
                aux2, 
                (const Float *)aux1, 
                Arows, Acolumns, Avalues, 
                auxR,
                inneriteration_threshold, 
                inneriteration_print_modulo
                , precon
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
            ConjugateGradientSolverCSR_DiagonalPreconditioner( 
                L, 
                aux2, 
                (const Float *)aux1, 
                Arows, Acolumns, Avalues, 
                auxR,
                inneriteration_threshold,
                inneriteration_print_modulo,
                precon
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
        
        /* Print information */
        
        if( print_modulo > 0 and k % print_modulo == 0 ) 
            printf( "Hodge Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
                    k, N, (long double)(Md_r), (long double) threshold*threshold );
        
        /* Check whether res is small */
                
        bool residualenergy_is_small = absolute(Md_r) < threshold*threshold;
        
        if( residualenergy_is_small )
            break;
        
        
        
        bool denominator_is_unreasonable = not std::isfinite(Md_Md) or Md_Md < 0.;
        bool denominator_is_small    = sqrt(absolute(Md_Md)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) {
            printf( "Gradient double energy is unreasonable with %.9Le\n", (long double)Md_Md );
            break;
        }
        
        if( denominator_is_small ) {
            printf( "Gradient double energy is small with %.9Le while precon-residual is %.9Le vs %.9Le\n", 
                    (long double)Md_Md, (long double)Md_r, (long double)threshold*threshold );
            break;
        }


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
        ConjugateGradientSolverCSR_DiagonalPreconditioner( 
            L, 
            aux2, 
            (const Float *)aux1, 
            Arows, Acolumns, Avalues, 
            auxR,
            inneriteration_threshold,
            inneriteration_print_modulo
            , precon
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
        
        k++;
        
    }
    
    if( print_modulo >= 0 ) 
        printf("Hodge Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
               k, N, (long double)(Md_r), (long double) threshold*threshold );

    
    free(  dir );
    free( Mdir );
    free( Mres );

    free( aux1 );
    free( aux2 );
    free( auxR );
    
    free( vil );
    
    free( precon );

}
 































void HodgeConjugateResidualSolverCSR_SSOR( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    Float* res,
    Float threshold,
    int print_modulo,
    Float inneriteration_threshold,
    int inneriteration_print_modulo
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
    assert( threshold > 0 );
//     assert( print_modulo >= 0 );
    
    Float* __restrict__  dir = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ Mdir = (Float*)malloc( sizeof(Float) * N );
    Float* __restrict__ Mres = (Float*)malloc( sizeof(Float) * N );
    
    Float* __restrict__ aux1 = (Float*)malloc( sizeof(Float) * L );
    Float* __restrict__ aux2 = (Float*)malloc( sizeof(Float) * L );
    Float* __restrict__ auxR = (Float*)malloc( sizeof(Float) * L );
    
    Float* __restrict__  vil = (Float*)malloc( sizeof(Float) * N );
    
    Float* __restrict__  diagonal = (Float*)malloc( sizeof(Float) * L );
    
    assert(  dir );
    assert( Mdir );
    assert( Mres );
    
    assert( aux1 );
    assert( aux2 );
    assert( auxR );
    
    assert( vil );
    
    assert( diagonal );
    
    
    #pragma omp parallel for
    for( int c = 0; c < L; c++ ) {
        
        diagonal[c] = 0.;
        
        for( int d = Arows[c]; d < Arows[c+1]; d++ )
            if( Acolumns[d] == c ) 
                diagonal[c] += Avalues[ d ];
            
//         assert( diagonal[c] >= 0. );
        
    }
    
    
    
    Float Md_r  = notanumber;
    Float Md_Md = notanumber;
    
    int k = 0;
    
    while( k < N ){
        
        bool restart_condition = ( k == 0 ); // or k % 1000 == 0;
        
        bool residualenergy_seems_small = false;//absolute(Md_r) < threshold*threshold;

        if( restart_condition or residualenergy_seems_small ) {
            
            #pragma omp parallel for
            for( int c = 0; c < L; c++ ) {
                
                aux1[c] = 0.;
                
                for( int d = Btrows[c]; d < Btrows[c+1]; d++ )
                    aux1[c] += Btvalues[ d ] * x[ Btcolumns[d] ];
                
            }
            
            for( int c = 0; c < L; c++ ) aux2[c] = 0.;
            ConjugateGradientSolverCSR_SSOR( 
                L, 
                aux2, 
                (const Float *)aux1, 
                Arows, Acolumns, Avalues, 
                auxR,
                inneriteration_threshold,
                inneriteration_print_modulo
                , diagonal, 1.
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
            ConjugateGradientSolverCSR_SSOR( 
                L, 
                aux2, 
                (const Float *)aux1, 
                Arows, Acolumns, Avalues, 
                auxR,
                inneriteration_threshold,
                inneriteration_print_modulo,
                diagonal, 1.
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
        
        /* Print information */
        
        if( print_modulo > 0 and k % print_modulo == 0 ) 
            printf("Hodge Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
                   k, N, (long double)(Md_r), (long double) threshold*threshold );
        
        /* Check whether res is small */
                
        bool residualenergy_is_small = absolute(Md_r) < threshold*threshold;
        
        if( residualenergy_is_small )
            break;
        
        
        
        bool denominator_is_unreasonable = not std::isfinite(Md_Md) or Md_Md < 0.;
        bool denominator_is_small    = sqrt(absolute(Md_Md)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) {
            printf( "Gradient double energy is unreasonable with %.9Le\n", (long double)Md_Md );
            break;
        }
        
        if( denominator_is_small ) {
            printf( "Gradient double energy is small with %.9Le while precon-residual is %.9Le vs %.9Le\n", 
                    (long double)Md_Md, (long double)Md_r, (long double)threshold*threshold );
            break;
        }


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
        ConjugateGradientSolverCSR_SSOR( 
            L, 
            aux2, 
            (const Float *)aux1, 
            Arows, Acolumns, Avalues, 
            auxR,
            inneriteration_threshold,
            inneriteration_print_modulo
            , diagonal, 1.0
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
        
        
        k++;
        
    }
    
    if( print_modulo >= 0 ) 
        printf("Hodge Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
               k, N, (long double)(Md_r), (long double) threshold*threshold );

    
    free(  dir );
    free( Mdir );
    free( Mres );

    free( aux1 );
    free( aux2 );
    free( auxR );
    
    free( vil );
    
    free( diagonal );

}
 
































void HodgeConjugateResidualSolverCSR_textbook( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    Float* res,
    Float threshold,
    int print_modulo,
    Float inneriteration_threshold,
    int inneriteration_print_modulo
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
    assert( threshold > 0 );
//     assert( print_modulo >= 0 );
    
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
    
    
    Float Mr_r  = notanumber;
    Float Md_Md = notanumber;
    
    int k = 0;
    
    while( k < N ){
        
        bool restart_condition = ( k == 0 ); // or k % 1000 == 0;
        
        bool residualenergy_seems_small = false;//absolute(Mr_r) < threshold*threshold;

        if( restart_condition or residualenergy_seems_small ) {
            
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
                inneriteration_threshold,
                inneriteration_print_modulo
            );
            
            #pragma omp parallel for
            for( int c = 0; c < N; c++ ) {
                
                res[c] = b[c];
                
                for( int d = Brows[c]; d < Brows[c+1]; d++ )
                    res[c] -= Bvalues[ d ] * aux2[ Bcolumns[d] ];
                
                dir[c] = res[c]; 
            
            }
            
            Mr_r  = 0.;
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
                inneriteration_threshold,
                inneriteration_print_modulo
            );
            
            #pragma omp parallel for reduction(+:Mr_r,Md_Md)
            for( int c = 0; c < N; c++ ) {
                
                Mdir[c] = 0.;
                
                for( int d = Brows[c]; d < Brows[c+1]; d++ )
                    Mdir[c] += Bvalues[ d ] * aux2[ Bcolumns[d] ];
                
                Mres[c] = Mdir[c];
                
                Mr_r  += Mres[c] *  res[c];
                Md_Md += Mdir[c] * Mdir[c];
                
            }
            
            
        }
        
        /* Print information */
        
        if( print_modulo > 0 and k % print_modulo == 0 ) 
            printf("Hodge Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
                   k, N, (long double)(Mr_r), (long double) threshold*threshold );
        
        /* Check whether res is small */
                
        bool residualenergy_is_small = absolute(Mr_r) < threshold*threshold;
        
        if( residualenergy_is_small )
            break;
        
        
        
        bool denominator_is_unreasonable = not std::isfinite(Md_Md) or Md_Md < 0.;
        bool denominator_is_small    = sqrt(absolute(Md_Md)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) {
            if( print_modulo >= 0 ) printf( "Gradient double energy is unreasonable with %.9Le\n", (long double)Md_Md );
            break;
        }
        
        if( denominator_is_small ) {
            if( print_modulo >= 0 ) printf( "Gradient double energy is small with %.9Le while precon-residual is %.9Le vs %.9Le\n", 
                                    (long double)Md_Md, (long double)Mr_r, (long double)threshold*threshold );
            break;
        }


        /* now the main work of the entire algorithm */
        
        Float alpha = Mr_r / Md_Md;
        
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
            inneriteration_threshold,
            inneriteration_print_modulo
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
        
        Float beta = new_Mr_r / Mr_r;
        
        Md_Md = 0.;
        #pragma omp parallel for reduction(+:Md_Md)
        for( int c = 0; c < N; c++ )
        {
            
             dir[c] =  res[c] + beta *  dir[c];
            
            Mdir[c] = Mres[c] + beta * Mdir[c];
            
            Md_Md += Mdir[c] * Mdir[c];
        }
        
        Mr_r = new_Mr_r;
        
        
        k++;
        
    }
    
    if( print_modulo >= 0 ) 
        printf("Hodge Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
               k, N, (long double)(Mr_r), (long double) threshold*threshold );

    
    free(  dir );
    free( Mdir );
    free( Mres );

    free( aux1 );
    free( aux2 );
    free( auxR );
    
    free( vil );

}




 

 
#endif
