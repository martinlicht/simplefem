
#include "systemsparsesolver.hpp"









void HodgeConjugateResidualSolverCSR( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    const int* __restrict__  Crows, const int* __restrict__  Ccolumns, const Float* __restrict__  Cvalues, 
    Float* res,
    Float threshold,
    int print_modulo,
    Float inneriteration_threshold,
    int inneriteration_print_modulo
) {
    HodgeConjugateResidualSolverCSR_SSOR( 
        N, 
        L, 
        x, 
        b, 
        Arows,   Acolumns,  Avalues, 
        Brows,   Bcolumns,  Bvalues, 
        Btrows, Btcolumns, Btvalues, 
        Crows,   Ccolumns,  Cvalues, 
        res,
        threshold,
        print_modulo,
        inneriteration_threshold,
        inneriteration_print_modulo
    );
}
 

void HodgeConjugateResidualSolverCSR_diagonal( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    const int* __restrict__  Crows, const int* __restrict__  Ccolumns, const Float* __restrict__  Cvalues, 
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
    
    Float* __restrict__  dir = new (std::nothrow) Float[N];
    Float* __restrict__ Mdir = new (std::nothrow) Float[N];
    Float* __restrict__ Mres = new (std::nothrow) Float[N];
    
    Float* __restrict__ aux1 = new (std::nothrow) Float[L];
    Float* __restrict__ aux2 = new (std::nothrow) Float[L];
    Float* __restrict__ auxR = new (std::nothrow) Float[L];
    
    Float* __restrict__  vil = new (std::nothrow) Float[N];
    
    Float* __restrict__  precon = new (std::nothrow) Float[L];
    
    assert(  dir );
    assert( Mdir );
    assert( Mres );
    
    assert( aux1 );
    assert( aux2 );
    assert( auxR );
    
    assert( vil );
    
    assert( precon );
    
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif 
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
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif 
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
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif
            for( int c = 0; c < N; c++ ) {
                
                res[c] = b[c];
                
                for( int d = Brows[c]; d < Brows[c+1]; d++ )
                    res[c] -= expected_sign_of_A * Bvalues[ d ] * aux2[ Bcolumns[d] ];
                
                for( int d = Crows[c]; d < Crows[c+1]; d++ )
                    res[c] -= Cvalues[ d ] *    x[ Ccolumns[d] ];
                
                dir[c] = res[c]; 
            
            }
            
            Md_r  = 0.;
            Md_Md = 0.;
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif 
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
            
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(+:Md_r,Md_Md)
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                Mdir[c] = 0.;
                
                for( int d = Brows[c]; d < Brows[c+1]; d++ )
                    Mdir[c] += expected_sign_of_A * Bvalues[ d ] * aux2[ Bcolumns[d] ];
                
                for( int d = Crows[c]; d < Crows[c+1]; d++ )
                    Mdir[c] += Cvalues[ d ] *  dir[ Ccolumns[d] ];
                
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
        
        #if defined(_OPENMP)
        #pragma omp parallel for
        #endif 
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
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:new_Mr_r)
        // #pragma omp parallel for 
        #endif
        for( int c = 0; c < N; c++ ) {
            
            vil[c] = 0.;
            
            for( int d = Brows[c]; d < Brows[c+1]; d++ )
                vil[c] += expected_sign_of_A * Bvalues[ d ] * aux2[ Bcolumns[d] ];
            
            for( int d = Crows[c]; d < Crows[c+1]; d++ )
                vil[c] += Cvalues[ d ] * Mdir[ Ccolumns[d] ];
            
        // }
                    
        // #if defined(_OPENMP)
        // #pragma omp parallel for reduction(+:new_Ar_r)
        // #endif 
        // for( int c = 0; c < N; c++ )
        // {
        
            x[c]    =   x[c] + alpha *  dir[c];
            
            res[c]  = res[c] - alpha * Mdir[c];
            
            Mres[c] = Mres[c] - alpha * vil[c];
            
            new_Mr_r += Mres[c] * res[c];
            
        }
        
        Float beta = new_Mr_r / Md_r;
        
        Md_Md = 0.;
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:Md_Md)
        #endif 
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
        printf("Final Hodge Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
               k, N, (long double)(Md_r), (long double) threshold*threshold );

    
    delete[] (  dir );
    delete[] ( Mdir );
    delete[] ( Mres );

    delete[] ( aux1 );
    delete[] ( aux2 );
    delete[] ( auxR );
    
    delete[] ( vil );
    
    delete[] ( precon );

}
 































void HodgeConjugateResidualSolverCSR_SSOR( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    const int* __restrict__  Crows, const int* __restrict__  Ccolumns, const Float* __restrict__  Cvalues, 
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
    
    Float* __restrict__  dir = new (std::nothrow) Float[N];
    Float* __restrict__ Mdir = new (std::nothrow) Float[N];
    Float* __restrict__ Mres = new (std::nothrow) Float[N];
    
    Float* __restrict__ aux1 = new (std::nothrow) Float[L];
    Float* __restrict__ aux2 = new (std::nothrow) Float[L];
    Float* __restrict__ auxR = new (std::nothrow) Float[L];
    
    Float* __restrict__  vil = new (std::nothrow) Float[N];
    
    Float* __restrict__  diagonal = new (std::nothrow) Float[L];
    
    assert(  dir );
    assert( Mdir );
    assert( Mres );
    
    assert( aux1 );
    assert( aux2 );
    assert( auxR );
    
    assert( vil );
    
    assert( diagonal );
    
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif 
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
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif 
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
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                res[c] = b[c];
                
                for( int d = Brows[c]; d < Brows[c+1]; d++ )
                    res[c] -= expected_sign_of_A * Bvalues[ d ] * aux2[ Bcolumns[d] ];
                
                for( int d = Crows[c]; d < Crows[c+1]; d++ )
                    res[c] -= Cvalues[ d ] *   x[ Ccolumns[d] ];
                
                dir[c] = res[c]; 
            
            }
            
            Md_r  = 0.;
            Md_Md = 0.;
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif
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
            
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(+:Md_r,Md_Md)
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                Mdir[c] = 0.;
                
                for( int d = Brows[c]; d < Brows[c+1]; d++ )
                    Mdir[c] += expected_sign_of_A * Bvalues[ d ] * aux2[ Bcolumns[d] ];
                
                for( int d = Crows[c]; d < Crows[c+1]; d++ )
                    Mdir[c] += Cvalues[ d ] *  dir[ Ccolumns[d] ];
                
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
        
        #if defined(_OPENMP)
        #pragma omp parallel for
        #endif 
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
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:new_Mr_r)
        // #pragma omp parallel for 
        #endif 
        for( int c = 0; c < N; c++ ) {
            
            vil[c] = 0.;
            
            for( int d = Brows[c]; d < Brows[c+1]; d++ )
                vil[c] += expected_sign_of_A * Bvalues[ d ] * aux2[ Bcolumns[d] ];
            
            for( int d = Crows[c]; d < Crows[c+1]; d++ )
                vil[c] += Cvalues[ d ] * Mdir[ Ccolumns[d] ];
            
        // }
                    
        // #if defined(_OPENMP)
        // #pragma omp parallel for reduction(+:new_Ar_r)
        // #endif 
        // for( int c = 0; c < N; c++ )
        // {
        
            x[c]    =   x[c] + alpha *  dir[c];
            
            res[c]  = res[c] - alpha * Mdir[c];
            
            Mres[c] = Mres[c] - alpha * vil[c];
            
            new_Mr_r += Mres[c] * res[c];
            
        }
        
        Float beta = new_Mr_r / Md_r;
        
        Md_Md = 0.;
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:Md_Md)
        #endif 
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
        printf("Final Hodge Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
               k, N, (long double)(Md_r), (long double) threshold*threshold );

    
    delete[] (  dir );
    delete[] ( Mdir );
    delete[] ( Mres );

    delete[] ( aux1 );
    delete[] ( aux2 );
    delete[] ( auxR );
    
    delete[] ( vil );
    
    delete[] ( diagonal );

}
 
































void HodgeConjugateResidualSolverCSR_textbook( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    const int* __restrict__  Crows, const int* __restrict__  Ccolumns, const Float* __restrict__  Cvalues, 
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
    
    Float* __restrict__  dir = new (std::nothrow) Float[N];
    Float* __restrict__ Mdir = new (std::nothrow) Float[N];
    Float* __restrict__ Mres = new (std::nothrow) Float[N];
    
    Float* __restrict__ aux1 = new (std::nothrow) Float[L];
    Float* __restrict__ aux2 = new (std::nothrow) Float[L];
    Float* __restrict__ auxR = new (std::nothrow) Float[L];
    
    Float* __restrict__  vil = new (std::nothrow) Float[N];
    
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
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif 
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
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                res[c] = b[c];
                
                for( int d = Brows[c]; d < Brows[c+1]; d++ ) 
                    res[c] -= expected_sign_of_A * Bvalues[ d ] * aux2[ Bcolumns[d] ];
                
                for( int d = Crows[c]; d < Crows[c+1]; d++ ) 
                    res[c] -= Cvalues[ d ] *    x[ Ccolumns[d] ];
                
                dir[c] = res[c]; 
            
            }
            
            Mr_r  = 0.;
            Md_Md = 0.;
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif 
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
            
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(+:Mr_r,Md_Md)
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                Mdir[c] = 0.;
                
                for( int d = Brows[c]; d < Brows[c+1]; d++ )
                    Mdir[c] += expected_sign_of_A * Bvalues[ d ] * aux2[ Bcolumns[d] ];
                
                for( int d = Crows[c]; d < Crows[c+1]; d++ )
                    Mdir[c] += Cvalues[ d ] *  dir[ Ccolumns[d] ];
                
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
        
        #if defined(_OPENMP)
        #pragma omp parallel for
        #endif 
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
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:new_Mr_r)
        // #pragma omp parallel for 
        #endif 
        for( int c = 0; c < N; c++ ) {
            
            vil[c] = 0.;
            
            for( int d = Brows[c]; d < Brows[c+1]; d++ )
                vil[c] += expected_sign_of_A * Bvalues[ d ] * aux2[ Bcolumns[d] ];
            
            for( int d = Crows[c]; d < Crows[c+1]; d++ )
                vil[c] += Cvalues[ d ] * Mdir[ Ccolumns[d] ];
            
        // }
                    
        #if defined(_OPENMP)
        // #pragma omp parallel for reduction(+:new_Ar_r)
        #endif 
        // for( int c = 0; c < N; c++ )
        // {
        
            x[c]    =   x[c] + alpha *  dir[c];
            
            res[c]  = res[c] - alpha * Mdir[c];
            
            Mres[c] = Mres[c] - alpha * vil[c];
            
            new_Mr_r += Mres[c] * res[c];
            
        }
        
        Float beta = new_Mr_r / Mr_r;
        
        Md_Md = 0.;
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:Md_Md)
        #endif 
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
        printf("Final Hodge Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
               k, N, (long double)(Mr_r), (long double) threshold*threshold );

    
    delete[] (  dir );
    delete[] ( Mdir );
    delete[] ( Mres );

    delete[] ( aux1 );
    delete[] ( aux2 );
    delete[] ( auxR );
    
    delete[] ( vil );

}




 

