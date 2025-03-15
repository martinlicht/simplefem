
#include <cmath>

#include <new>
#include <utility>

#include "systemsparsesolver.hpp"

#include "sparsesolver.hpp"

static const bool csrsys_restart_on_full_dimension = false;
static const bool csrsys_restart_before_finish     = false;





int HodgeConjugateResidualSolverCSR( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    const int* __restrict__  Crows, const int* __restrict__  Ccolumns, const Float* __restrict__  Cvalues, 
    Float* res,
    Float precision,
    int print_modulo,
    Float inneriteration_precision,
    int inneriteration_print_modulo
) {
    return 
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
        precision,
        print_modulo,
        inneriteration_precision,
        inneriteration_print_modulo
    );
}
 

int HodgeConjugateResidualSolverCSR_diagonal( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    const int* __restrict__  Crows, const int* __restrict__  Ccolumns, const Float* __restrict__  Cvalues, 
    Float* res,
    Float precision,
    int print_modulo,
    Float inneriteration_precision,
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
    assert( precision > 0 );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = maximum( desired_precision, precision * sqrt(tolerance) );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    const bool do_print_breakdown = print_modulo >=  0;
    const bool do_print_warning   = print_modulo >=  0;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
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
    
    int K = 0;
    
    if( do_print_begin ) LOGPRINTF( "(%d/%d)     BEGIN: Hodge Conjugate Residual (diag) CSR\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csrsys_restart_on_full_dimension and K % N == 0 );
        
        bool residualenergy_seems_small = ( K != 0 ) and absolute(Md_r) < tolerance*tolerance;
        // bool residualenergy_seems_small = false;

        if( restart_condition or ( residualenergy_seems_small and csrsys_restart_before_finish ) ) UNLIKELY {
            
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
                inneriteration_precision, 
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
                inneriteration_precision,
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
            
            if( do_print_restart ) 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Md_r), (double)(safedouble) tolerance );

        }
        
        /* Print information */

        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );
        
        if( print_condition and do_print_interim ) UNLIKELY 
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Md_r), (double)(safedouble) tolerance );
        
        /* Check whether res is small */
        
        bool residualenergy_is_unreasonable = not std::isfinite(Md_r) or Md_r < 0.;
        bool residualenergy_is_small = absolute(Md_r) < tolerance*tolerance;
        
        if( residualenergy_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Residual energy is unreasonable with %.9le\n", K, N, (double)(safedouble)Md_r );
            break;
        }

        if( residualenergy_is_small ) UNLIKELY 
            break;
        
        
        
        bool denominator_is_unreasonable = not std::isfinite(Md_Md) or Md_Md < 0.;
        bool denominator_is_small        = std::sqrt(absolute(Md_Md)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient double energy is unreasonable with %.9le\n", K, N, (double)(safedouble)Md_Md );
            break;
        }
        
        if( denominator_is_small ) UNLIKELY {
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Md_r), (double)(safedouble) tolerance );
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient double energy is small with %.9le\n", K, N, (double)(safedouble)Md_Md );
            break;
        }


        /* now the main work of the entire algorithm */
        
        Float alpha = Md_r / Md_Md;
        
        Float new_Mr_r = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for
        #endif 
        for( int c = 0; c < L; c++ ) {
            
            aux1[c] = 0.;
            
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
            inneriteration_precision,
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
        
        K++;
        
    }
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(Md_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble) std::sqrt(Md_r), (double)(safedouble) tolerance );

    
    delete[] (  dir );
    delete[] ( Mdir );
    delete[] ( Mres );

    delete[] ( aux1 );
    delete[] ( aux2 );
    delete[] ( auxR );
    
    delete[] ( vil );
    
    delete[] ( precon );

    return K;

}
 































int HodgeConjugateResidualSolverCSR_SSOR( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    const int* __restrict__  Crows, const int* __restrict__  Ccolumns, const Float* __restrict__  Cvalues, 
    Float* res,
    Float precision,
    int print_modulo,
    Float inneriteration_precision,
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
    assert( precision > 0 );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = maximum( desired_precision, precision * sqrt(tolerance) );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    // const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    const bool do_print_breakdown = print_modulo >=  0;
    // const bool do_print_warning   = print_modulo >=  0;
    // const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
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
    
    int K = 0;
    
    if( do_print_begin ) LOGPRINTF( "(%d/%d)     BEGIN: Hodge Conjugate Residual (SSOR) CSR\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csrsys_restart_on_full_dimension and K % N == 0 );
        
        bool residualenergy_seems_small = ( K != 0 ) and absolute(Md_r) < tolerance*tolerance;
        // bool residualenergy_seems_small = false;

        if( restart_condition or ( residualenergy_seems_small and csrsys_restart_before_finish ) ) UNLIKELY {
            
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
                inneriteration_precision,
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
                
                aux1[c] = 0.;
                
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
                inneriteration_precision,
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
            
            if( do_print_restart ) 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Md_r), (double)(safedouble) tolerance );

        }
        
        /* Print information */
        
        if( print_modulo > 0 and K % print_modulo == 0 ) UNLIKELY 
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Md_r), (double)(safedouble) tolerance );
        
        /* Check whether res is small */
                
        bool residualenergy_is_unreasonable = not std::isfinite(Md_r) or Md_r < 0.;

        if( residualenergy_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Residual energy is unreasonable with %.9le\n", K, N, (double)(safedouble)Md_r );
            break;
        }

        bool residualenergy_is_small = absolute(Md_r) < tolerance*tolerance;
        
        if( residualenergy_is_small ) UNLIKELY 
            break;
        
        
        
        bool denominator_is_unreasonable = not std::isfinite(Md_Md) or Md_Md < 0.;
        bool denominator_is_small        = std::sqrt(absolute(Md_Md)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) UNLIKELY {
            LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient double energy is unreasonable with %.9le\n", K, N, (double)(safedouble)Md_Md );
            break;
        }
        
        if( denominator_is_small ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Md_r), (double)(safedouble) tolerance );
            LOGPRINTF( "(%d/%d)   WARNING: Gradient double energy is small with %.9le\n", K, N, (double)(safedouble)Md_Md );
            break;
        }


        /* now the main work of the entire algorithm */
        
        Float alpha = Md_r / Md_Md;
        
        Float new_Mr_r = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for
        #endif 
        for( int c = 0; c < L; c++ ) {
            
            aux1[c] = 0.;
            
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
            inneriteration_precision,
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
        
        
        K++;
        
    }
    
    if( print_modulo >= 0 ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(Md_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble) std::sqrt(Md_r), (double)(safedouble) tolerance );

    
    delete[] (  dir );
    delete[] ( Mdir );
    delete[] ( Mres );

    delete[] ( aux1 );
    delete[] ( aux2 );
    delete[] ( auxR );
    
    delete[] ( vil );
    
    delete[] ( diagonal );

    return K;

}
 
































int HodgeConjugateResidualSolverCSR_textbook( 
    const int N, 
    const int L, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    const int* __restrict__  Crows, const int* __restrict__  Ccolumns, const Float* __restrict__  Cvalues, 
    Float* res,
    Float precision,
    int print_modulo,
    Float inneriteration_precision,
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
    assert( precision > 0 );

    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = maximum( desired_precision, precision * sqrt(tolerance) );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    const bool do_print_breakdown = print_modulo >=  0;
    const bool do_print_warning   = print_modulo >=  0;
    const bool do_print_finish    = print_modulo >= -1;

    /* Build up data */

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
    
    int K = 0;
    
    if( do_print_begin ) LOGPRINTF( "(%d/%d)     BEGIN: Hodge Conjugate Residual (textbook) CSR\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csrsys_restart_on_full_dimension and K % N == 0 );
        
        bool residualenergy_seems_small = ( K != 0 ) and absolute(Mr_r) < tolerance*tolerance;

        if( restart_condition or ( residualenergy_seems_small and csrsys_restart_before_finish ) ) {
            
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
                inneriteration_precision,
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
                
                aux1[c] = 0.;
                
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
                inneriteration_precision,
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

            if( do_print_restart ) 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Mr_r), (double)(safedouble) tolerance );
            
        }
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );
        
        if( print_condition and do_print_interim ) 
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Mr_r), (double)(safedouble) tolerance );
        
        /* Check whether res is small */
                
        bool residualenergy_is_unreasonable = not std::isfinite(Mr_r) or Mr_r < 0.;
        
        if( residualenergy_is_unreasonable ) {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Residual energy is unreasonable with %.9le\n", K, N, (double)(safedouble)Mr_r );
            break;
        }

        bool residualenergy_is_small = absolute(Mr_r) < tolerance*tolerance;
        
        if( residualenergy_is_small )
            break;
        
        
        
        bool denominator_is_unreasonable = not std::isfinite(Md_Md) or Md_Md < 0.;
        bool denominator_is_small        = std::sqrt(absolute(Md_Md)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient double energy is unreasonable with %.9le\n", K, N, (double)(safedouble)Md_Md );
            break;
        }
        
        if( denominator_is_small ) {
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Mr_r), (double)(safedouble) tolerance );
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient double energy is small with %.9le\n", K, N, (double)(safedouble)Md_Md );
            break;
        }


        /* now the main work of the entire algorithm */
        
        Float alpha = Mr_r / Md_Md;
        
        Float new_Mr_r = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for
        #endif 
        for( int c = 0; c < L; c++ ) {
            
            aux1[c] = 0.;
            
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
            inneriteration_precision,
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
        
        
        K++;
        
    }
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(Mr_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble) std::sqrt(Mr_r), (double)(safedouble) tolerance );

    
    delete[] (  dir );
    delete[] ( Mdir );
    delete[] ( Mres );

    delete[] ( aux1 );
    delete[] ( aux2 );
    delete[] ( auxR );
    
    delete[] ( vil );

    return K;

}




 





int HodgeHerzogSoodhalterMethod( 
    const int dimension_A, 
    const int dimension_C, 
    Float* __restrict__ x_A, 
    Float* __restrict__ x_C, 
    const Float* __restrict__ b_A, 
    const Float* __restrict__ b_C, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    const int* __restrict__  Crows, const int* __restrict__  Ccolumns, const Float* __restrict__  Cvalues, 
    Float precision,
    int print_modulo,
    const int* __restrict__ PArows, const int* __restrict__ PAcolumns, const Float* __restrict__ PAvalues, 
    const int* __restrict__ PCrows, const int* __restrict__ PCcolumns, const Float* __restrict__ PCvalues, 
    Float inneriteration_precision,
    int inneriteration_print_modulo
) {

    assert( dimension_A > 0 );
    assert( dimension_C > 0 );
    assert( x_A );
    assert( x_C );
    assert( b_A );
    assert( b_C );
    assert(  Arows );
    assert(  Acolumns );
    assert(  Avalues );
    assert(  Brows );
    assert(  Bcolumns );
    assert(  Bvalues );
    assert( Btrows );
    assert( Btcolumns );
    assert( Btvalues );
    assert( precision > 0 );
    // assert( PArows ); 
    // assert( PAcolumns );
    // assert( PAvalues );
    // assert( PCrows );
    // assert( PCcolumns );
    // assert( PCvalues );
    assert( inneriteration_precision > 0 );
//     assert( inneriteration_print_modulo >= 0 );
    
    Float tolerance = 0.;
    for( int i = 0; i < dimension_A; i++ ) tolerance += b_A[i]*b_A[i];
    for( int i = 0; i < dimension_C; i++ ) tolerance += b_C[i]*b_C[i];
    tolerance = maximum( desired_precision, precision * sqrt(tolerance) );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    // const bool do_print_breakdown = print_modulo >=  0;
    // const bool do_print_warning   = print_modulo >=  0;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    Float* __restrict__ v0_A = new (std::nothrow) Float[ dimension_A ];
    Float* __restrict__ v1_A = new (std::nothrow) Float[ dimension_A ];
    Float* __restrict__ w0_A = new (std::nothrow) Float[ dimension_A ];
    Float* __restrict__ w1_A = new (std::nothrow) Float[ dimension_A ];
    Float* __restrict__  z_A = new (std::nothrow) Float[ dimension_A ];
    
    Float* __restrict__ v0_C = new (std::nothrow) Float[ dimension_C ];
    Float* __restrict__ v1_C = new (std::nothrow) Float[ dimension_C ];
    Float* __restrict__ w0_C = new (std::nothrow) Float[ dimension_C ];
    Float* __restrict__ w1_C = new (std::nothrow) Float[ dimension_C ];
    Float* __restrict__  z_C = new (std::nothrow) Float[ dimension_C ];
    
    Float* __restrict__ vn_A = new (std::nothrow) Float[ dimension_A ];
    Float* __restrict__ wn_A = new (std::nothrow) Float[ dimension_A ];
    Float* __restrict__ zn_A = new (std::nothrow) Float[ dimension_A ];
    
    Float* __restrict__ vn_C = new (std::nothrow) Float[ dimension_C ];
    Float* __restrict__ wn_C = new (std::nothrow) Float[ dimension_C ];
    Float* __restrict__ zn_C = new (std::nothrow) Float[ dimension_C ];
    
    Float* __restrict__  m_A = new (std::nothrow) Float[ dimension_A ];
    Float* __restrict__  m_C = new (std::nothrow) Float[ dimension_C ];
    
    Float* __restrict__  p_A = new (std::nothrow) Float[ dimension_A ];
    Float* __restrict__  p_C = new (std::nothrow) Float[ dimension_C ];
    
    
    Float mu_A = notanumber;
    Float mu_C = notanumber;
    
    Float gamma = notanumber;
    
    Float eta   = notanumber;
    Float eta_A = notanumber;
    Float eta_C = notanumber;
    
    Float s0 = notanumber;
    Float s1 = notanumber;
    Float c0 = notanumber;
    Float c1 = notanumber;
    
    int max_iteration_count = dimension_A + dimension_C;
    int recent_iteration_count = 0;


    bool precon_A_available = PArows and PAcolumns and PAvalues;
    bool precon_C_available = PCrows and PCcolumns and PCvalues;
    if( precon_A_available ) {
        assert( PArows ); 
        assert( PAcolumns );
        assert( PAvalues );
    }
    if( precon_C_available ) {
        assert( PCrows ); 
        assert( PCcolumns );
        assert( PCvalues );
    }
    
    if( do_print_begin ) LOGPRINTF( "(%d/%d)     BEGIN: Hodge Herzog-Soodhalter CSR\n", recent_iteration_count, max_iteration_count );
    
    if( precon_A_available and do_print_begin ) LOGPRINTF("      Preconditioner for A detected\n");
    if( precon_C_available and do_print_begin ) LOGPRINTF("      Preconditioner for C detected\n");

    while( recent_iteration_count < max_iteration_count ){
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( csrsys_restart_on_full_dimension and recent_iteration_count != 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and ( absolute(eta) < tolerance );
        
        if( restart_condition or ( residual_seems_small and csrsys_restart_before_finish ) ) {
            
            // 1
            for( int a = 0; a < dimension_A; a++ ) v0_A[a] = w0_A[a] = w1_A[a] = 0.;
            for( int c = 0; c < dimension_C; c++ ) v0_C[c] = w0_C[c] = w1_C[c] = 0.;
            
            // 2 
            for( int r = 0; r < dimension_A; r++ ) {
                v1_A[r] = b_A[r];
                for( int d = Arows[r]; d < Arows[r+1]; d++ )
                    v1_A[r] += expected_sign_of_A * Avalues[d] * x_A[ Acolumns[d] ];
                for( int d = Btrows[r]; d < Btrows[r+1]; d++ )
                    v1_A[r] -= Btvalues[d] * x_C[ Btcolumns[d] ];
            }
            for( int r = 0; r < dimension_C; r++ ) {
                v1_C[r] = b_C[r];
                for( int d = Brows[r]; d < Brows[r+1]; d++ )
                    v1_C[r] -= Bvalues[d] * x_A[ Bcolumns[d] ];
                for( int d = Crows[r]; d < Crows[r+1]; d++ )
                    v1_C[r] -= Cvalues[d] * x_C[ Ccolumns[d] ];
            }

            // In case we don't use preconditioners 
            for( int r = 0; r < dimension_A; r++ ) z_A[r] = v1_A[r]; 
            for( int r = 0; r < dimension_C; r++ ) z_C[r] = v1_C[r]; 
            
            // if preconditioners are available, we simply ditch the write above and compute z_A and z_C again from v1_A and v1_C, respectively.
            // z_A = PAinv * v1_A;
            // z_C = PCinv * v1_C;
            
            if( precon_A_available )
            {
                ConjugateGradientSolverCSR( 
                    dimension_A, 
                    z_A, 
                    v1_A, 
                    PArows, PAcolumns, PAvalues,
                    wn_A, // we recycle this memory 
                    inneriteration_precision,
                    inneriteration_print_modulo
                );
            }

            if( precon_C_available )
            {
                ConjugateGradientSolverCSR( 
                    dimension_C, 
                    z_C, 
                    v1_C, 
                    PCrows, PCcolumns, PCvalues,
                    wn_C, // we recycle this memory 
                    inneriteration_precision,
                    inneriteration_print_modulo
                );
            }
            
            
            // 3 -- 6
            Float v1_z_A = 0.;
            Float v1_z_C = 0.;
            for( int a = 0; a < dimension_A; a++ ) v1_z_A += v1_A[a] * z_A[a];
            for( int c = 0; c < dimension_C; c++ ) v1_z_C += v1_C[c] * z_C[c];
            
            Assert( v1_z_A + v1_z_C >= 0., v1_z_A, v1_z_C );
            gamma = std::sqrt( v1_z_A + v1_z_C );
            Assert( gamma > 0., recent_iteration_count, gamma );
            
            Assert( gamma > 0., v1_z_A, v1_z_C );
            
            for( int a = 0; a < dimension_A; a++ ){
                v1_A[a] /= gamma; z_A[a] /= gamma; 
                m_A[a] = v1_A[a];
            }

            for( int c = 0; c < dimension_C; c++ ){
                v1_C[c] /= gamma; z_C[c] /= gamma; 
                m_C[c] = v1_C[c];
            }
            
            Float psi_A = v1_z_A / (gamma*gamma); 
            Float psi_C = v1_z_C / (gamma*gamma); 
            mu_A = psi_A; 
            mu_C = psi_C; 
            
            // 7 
            eta = gamma; 
            eta_A = gamma * std::sqrt( psi_A );
            eta_C = gamma * std::sqrt( psi_C );
            
            s0 = s1 = 0.;
            c0 = c1 = 1.;
            
            if( do_print_restart ) {
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) absolute(eta), (double)(safedouble)tolerance );
                LOGPRINTF( "(%d/%d)            Gamma: %.9le Eta_A %.9le Eta_C %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)eta_A, (double)(safedouble)eta_C, (double)(safedouble)gamma );
            }

        }
        
        bool residual_is_small = ( absolute(eta) < tolerance );
        
        if( residual_is_small )
            break;

            
        {
            
            // 9  
            Float z_p_A = 0.;
            Float z_p_C = 0.;
            
            for( int r = 0; r < dimension_A; r++ ) {
                p_A[r] = 0.;
                for( int d = Arows[r]; d < Arows[r+1]; d++ )
                    p_A[r] -= expected_sign_of_A * Avalues[d]  * z_A[ Acolumns[d] ];
                for( int d = Btrows[r]; d < Btrows[r+1]; d++ )
                    p_A[r] += Btvalues[d] * z_C[ Btcolumns[d] ];
                z_p_A += p_A[r] * z_A[r];
            }
            for( int r = 0; r < dimension_C; r++ ) {
                p_C[r] = 0.;
                for( int d = Brows[r]; d < Brows[r+1]; d++ )
                    p_C[r] += Bvalues[d] * z_A[ Bcolumns[d] ];
                for( int d = Crows[r]; d < Crows[r+1]; d++ )
                    p_C[r] += Cvalues[d] * z_C[ Ccolumns[d] ];
                z_p_C += p_C[r] * z_C[r];
            }
                
 
            Float delta = z_p_A + z_p_C;

            // 10 
            for( int a = 0; a < dimension_A; a++ ) vn_A[a] = p_A[a] - delta * v1_A[a] - gamma * v0_A[a];
            for( int c = 0; c < dimension_C; c++ ) vn_C[c] = p_C[c] - delta * v1_C[c] - gamma * v0_C[c];
            
            // for( int a = 0; a < dimension_A; a++ ) {
            //     LOG << "A " << a << " : " << v1_A[a] << nl;
            // }
            // for( int c = 0; c < dimension_C; c++ ) {
            //     LOG << "C " << c << " : " << v1_C[c] << nl;
            // }

            

            // 11
            // In case we don't use preconditioners 
            Float zn_vn_A = 0.;
            Float zn_vn_C = 0.;
            for( int a = 0; a < dimension_A; a++ ) {
                zn_A[a] = vn_A[a];
                zn_vn_A += zn_A[a] * vn_A[a];
            }
            for( int c = 0; c < dimension_C; c++ ) {
                zn_C[c] = vn_C[c];
                zn_vn_C += zn_C[c] * vn_C[c];
            }
            
            
            // if preconditioners are available, we simply ditch the write above and compute z_A and z_C again from v1_A and v1_C, respectively.
            // z_A = PAinv * v1_A;
            // z_C = PCinv * v1_C;
            
            if( precon_A_available )
            {
                ConjugateGradientSolverCSR( 
                    dimension_A, 
                    zn_A, 
                    vn_A, 
                    PArows, PAcolumns, PAvalues,
                    wn_A, // we recycle this memory 
                    inneriteration_precision,
                    inneriteration_print_modulo
                );

                zn_vn_A = 0.;
                for( int a = 0; a < dimension_A; a++ ) zn_vn_A += zn_A[a] * vn_A[a];
            }

            if( precon_C_available )
            {
                ConjugateGradientSolverCSR( 
                    dimension_C, 
                    zn_C, 
                    vn_C, 
                    PCrows, PCcolumns, PCvalues,
                    wn_C, // we recycle this memory 
                    inneriteration_precision,
                    inneriteration_print_modulo
                );

                zn_vn_C = 0.;
                for( int c = 0; c < dimension_C; c++ ) zn_vn_C += zn_C[c] * vn_C[c];
            }
            

            // 12 
            Assert( zn_vn_A + zn_vn_C >= 0., zn_vn_A + zn_vn_C );
            Float gamma_n = std::sqrt( zn_vn_A + zn_vn_C );
            Assert( gamma_n > 0., recent_iteration_count, gamma_n, zn_vn_A, zn_vn_C );
            
            
            // 13 -- 14 
            for( int a = 0; a < dimension_A; a++ ){ vn_A[a] /= gamma_n; zn_A[a] /= gamma_n; }
            for( int c = 0; c < dimension_C; c++ ){ vn_C[c] /= gamma_n; zn_C[c] /= gamma_n; }

            // 15 -- 18
            Float alpha_0 = c1 * delta - c0 * s1 * gamma;
            assert( alpha_0 * alpha_0 + gamma_n * gamma_n > 0. );
            Float alpha_1 = std::sqrt( alpha_0 * alpha_0 + gamma_n * gamma_n );
            Float alpha_2 = s1 * delta + c0 * c1 * gamma;
            Float alpha_3 = s0 * gamma;
 
            assert( alpha_1 > 0. );

            // 19
            Float cn = alpha_0 / alpha_1;
            Float sn = gamma_n / alpha_1;
            
            // 20 -- 21 
            Float theta_A = 0.;
            Float theta_C = 0.;
            Float psi_A = 0.;
            Float psi_C = 0.;
            for( int a = 0; a < dimension_A; a++ ){ theta_A += m_A[a] * zn_A[a]; psi_A += zn_A[a] * vn_A[a]; }
            for( int c = 0; c < dimension_C; c++ ){ theta_C += m_C[c] * zn_C[c]; psi_C += zn_C[c] * vn_C[c]; }
            
            // 22 -- 24
            for( int a = 0; a < dimension_A; a++ ){
                m_A [a] = - sn * m_A[a] + cn * vn_A[a];
                wn_A[a] = ( z_A[a] - alpha_3 * w0_A[a] - alpha_2 * w1_A[a] ) / alpha_1;
                x_A [a] = x_A[a] + cn * eta * wn_A[a];
            }

            for( int c = 0; c < dimension_C; c++ ){ 
                m_C [c] = - sn * m_C[c] + cn * vn_C[c];
                wn_C[c] = ( z_C[c] - alpha_3 * w0_C[c] - alpha_2 * w1_C[c] ) / alpha_1;
                x_C [c] = x_C[c] + cn * eta * wn_C[c];
            }
 
            // 25 -- 26
            mu_A = sn * sn * mu_A - 2 * sn * cn * theta_A + cn * cn * psi_A;
            mu_C = sn * sn * mu_C - 2 * sn * cn * theta_C + cn * cn * psi_C;

            // 27 -- 28
            eta   = -sn * eta;
            eta_A = eta * std::sqrt( psi_A );
            eta_C = eta * std::sqrt( psi_C );

            // swaps ... 
            std::swap( v0_A, v1_A ); std::swap( v1_A, vn_A );
            std::swap( v0_C, v1_C ); std::swap( v1_C, vn_C );
            std::swap( w0_A, w1_A ); std::swap( w1_A, wn_A );
            std::swap( w0_C, w1_C ); std::swap( w1_C, wn_C );
            
            std::swap( z_A, zn_A );
            std::swap( z_C, zn_C );
            
            gamma = gamma_n;
            
            c0 = c1; c1 = cn;
            s0 = s1; s1 = sn;

            bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );

            if( do_print_interim and print_condition ) {
                LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) absolute(eta), (double)(safedouble)tolerance );
                LOGPRINTF( "(%d/%d)            Gamma: %.9le Eta_A %.9le Eta_C %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)eta_A, (double)(safedouble)eta_C, (double)(safedouble)gamma );
            }

        }

        
        Float recent_deviation = eta;
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( do_print_interim and print_condition ) {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) recent_deviation, (double)(safedouble)tolerance );
            LOGPRINTF( "(%d/%d)            Gamma: %.9le Eta: %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)gamma, (double)(safedouble)eta );
        }
        
        recent_iteration_count++;
        
    }
    
    delete[] v0_A;
    delete[] v1_A;
    delete[] w0_A;
    delete[] w1_A;
    delete[]  z_A;
    
    delete[] v0_C;
    delete[] v1_C;
    delete[] w0_C;
    delete[] w1_C;
    delete[]  z_C;
    
    delete[] vn_A;
    delete[] wn_A;
    delete[] zn_A;
    
    delete[] vn_C;
    delete[] wn_C;
    delete[] zn_C;
    
    delete[]  m_A;
    delete[]  m_C;
    
    delete[]  p_A;
    delete[]  p_C;

    /* HOW DID WE FINISH ? */
    
    Float recent_deviation = absolute( eta );
        
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "   "Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, recent_deviation < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble)recent_deviation, (double)(safedouble)tolerance );

    return recent_iteration_count;
}
  


