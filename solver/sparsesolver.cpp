



#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <new>
#include <utility>

#include "sparsesolver.hpp"


#if defined(_OPENMP)
#include <omp.h>
#endif


static const bool csr_restart_on_full_dimension = false;
static const bool csr_restart_before_finish     = false;



int ConjugateGradientSolverCSR( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float precision,
    int print_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    assert( precision > 0 );
    // assert( print_modulo >= -1 );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = precision * sqrt( tolerance );
    tolerance = maximum( tolerance, desired_precision );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    const bool do_print_breakdown = print_modulo >= -2;
    const bool do_print_warning   = print_modulo >= -2;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    Float* __restrict__ direction = new (std::nothrow) Float[N]; std::fill_n(direction,N,notanumber);
    Float* __restrict__ auxiliary = new (std::nothrow) Float[N]; std::fill_n(auxiliary,N,notanumber);
    assert( direction );
    assert( auxiliary );
    
    Float r_r = notanumber;

    int K = 0;
    
    if( do_print_begin ) 
        LOGPRINTF( "(%d/%d)     BEGIN: Conjugate Gradient CSR\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csr_restart_on_full_dimension and K % N == 0 );
        
        bool residual_seems_small = ( K != 0 ) and absolute(r_r) < tolerance*tolerance;

        if( restart_condition or ( residual_seems_small and csr_restart_before_finish ) ) UNLIKELY {
            
            r_r = 0.;
        
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(+:r_r)
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                residual[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                direction[c] = residual[c]; // this line seems to slow down performance ....
                
                r_r += residual[c] * residual[c];

            }

            if( do_print_restart ) 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(r_r), (double)(safedouble) tolerance );

        
        }
        
        /* printing information */

        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );

        if( print_condition and do_print_interim ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(r_r), (double)(safedouble) tolerance );
        }

        /* Check whether residual is small */
                
        bool residual_is_small = absolute(r_r) < tolerance*tolerance;
        
        if( residual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        // NOTE The calculation of d_r is reduced to r_r, which is already known.
        
        Float d_Ad = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:d_Ad) //d_r,
        #endif 
        for( int c = 0; c < N; c++ )
        {
            auxiliary[c] = 0.;
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                auxiliary[c] += csrvalues[ d ] * direction[ csrcolumns[d] ];
                    
            d_Ad += direction[c] * auxiliary[c];
        }
        
        
        bool numerator_is_unreasonable   = not std::isfinite(r_r) or r_r < 0.;
        bool denominator_is_unreasonable = not std::isfinite(d_Ad) or d_Ad < 0.;
        bool denominator_is_small        = std::sqrt(absolute(d_Ad)) < machine_epsilon;
        
        if( numerator_is_unreasonable or denominator_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient energy is unreasonable with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        if( denominator_is_small ) UNLIKELY {
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(r_r), (double)(safedouble) tolerance );
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient energy is small with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        
        Float alpha = r_r / d_Ad;
    
        Float r_r_new = 0.;
        
        #if defined(_OPENMP)    
        #pragma omp parallel for reduction(+:r_r_new) //r_r_old
        #endif 
        for( int c = 0; c < N; c++ )
        {
            
            x[c] += alpha * direction[c];
            
            residual[c] -= alpha * auxiliary[c];
            
            r_r_new += residual[c] * residual[c];
        }
        
        Float beta = r_r_new / r_r;
        
        r_r = r_r_new;
        
        #if defined(_OPENMP)
        #pragma omp parallel 
        #endif 
        for( int c = 0; c < N; c++ )
            direction[c] = residual[c] + beta * direction[c];
        
        
        
        K++;
        
    }
    
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(r_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble) std::sqrt(r_r), (double)(safedouble) tolerance );

    
    delete[] ( direction ); 
    delete[] ( auxiliary );

    return K;

}



















int ConjugateGradientSolverCSR_DiagonalPreconditioner( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float precision,
    int print_modulo,
    const Float* __restrict__ precon
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    assert( precision > 0 );
    // assert( print_modulo >= -1 );
    assert( precon );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = precision * sqrt( tolerance );
    tolerance = maximum( tolerance, desired_precision );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    const bool do_print_breakdown = print_modulo >= -2;
    const bool do_print_warning   = print_modulo >= -2;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    Float* __restrict__ direction = new (std::nothrow) Float[N]; std::fill_n(direction,N,notanumber);
    Float* __restrict__ zirconium = new (std::nothrow) Float[N]; std::fill_n(zirconium,N,notanumber);
    Float* __restrict__ auxiliary = new (std::nothrow) Float[N]; std::fill_n(auxiliary,N,notanumber);
    assert( direction );
    assert( zirconium );
    assert( auxiliary );
    
    Float z_r = notanumber;

    int K = 0;
    
    if( do_print_begin ) 
        LOGPRINTF( "(%d/%d)     BEGIN: Conjugate Gradient CSR (Diag)\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csr_restart_on_full_dimension and K % N == 0 );
        
        bool preconresidual_seems_small = ( K != 0 ) and absolute(z_r) < tolerance*tolerance;

        if( restart_condition or ( preconresidual_seems_small and csr_restart_before_finish ) ) UNLIKELY {
            
            z_r = 0.;
        
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(+:z_r)
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                if( precon[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
                residual[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                zirconium[c] = precon[c] * residual[c];
                
//                 assert( precon[c] == 1. );
                
                direction[c] = zirconium[c];
                
                z_r += zirconium[c] * residual[c];
            
            }

            if( do_print_restart ) 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );

            
        }
        
        /* printing information */

        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );

        if( print_condition and do_print_interim ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );
        }

        /* Check whether residual is small */
                
        bool preconresidual_is_small = absolute(z_r) < tolerance*tolerance;
        
        if( preconresidual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        // NOTE The calculation of d_r is reduced to r_r, which is already known.
        
        Float d_Ad = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:d_Ad) //d_r,
        #endif 
        for( int c = 0; c < N; c++ )
        {
            if( precon[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            auxiliary[c] = 0.;
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                auxiliary[c] += csrvalues[ d ] * direction[ csrcolumns[d] ];
                    
            d_Ad += direction[c] * auxiliary[c];
        }
        
        
        bool numerator_is_unreasonable   = not std::isfinite(z_r) or z_r < 0.;
        bool denominator_is_unreasonable = not std::isfinite(d_Ad) or d_Ad < 0.;
        bool denominator_is_small        = std::sqrt(absolute(d_Ad)) < machine_epsilon;
        
        if( numerator_is_unreasonable or denominator_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient energy is unreasonable with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        if( denominator_is_small ) UNLIKELY {
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient energy is small with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        
        Float alpha = z_r / d_Ad;
    
        Float z_r_new = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:z_r_new) //r_r_old
        #endif 
        for( int c = 0; c < N; c++ )
        {
            
            if( precon[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            x[c] += alpha * direction[c];
            
            residual[c] -= alpha * auxiliary[c];
            
            zirconium[c] = precon[c] * residual[c];
            
            z_r_new += zirconium[c] * residual[c];
        }
        
        Float beta = z_r_new / z_r;
        
        z_r = z_r_new;
        
        #if defined(_OPENMP)
        #pragma omp parallel 
        #endif
        for( int c = 0; c < N; c++ ) {
            
            if( precon[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            direction[c] = zirconium[c] + beta * direction[c];
            
        }
        
        K++;
        
    }
    
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(z_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );

    
    delete[] ( direction ); 
    delete[] ( zirconium ); 
    delete[] ( auxiliary );

    return K;

}





























int ConjugateGradientSolverCSR_SSOR( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float precision,
    int print_modulo,
    const Float* __restrict__ diagonal,
    Float omega
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    assert( precision > 0 );
    // assert( print_modulo >= -1 );
    assert( diagonal );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = precision * sqrt( tolerance );
    tolerance = maximum( tolerance, desired_precision );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    const bool do_print_breakdown = print_modulo >= -2;
    const bool do_print_warning   = print_modulo >= -2;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    Float* __restrict__ direction = new (std::nothrow) Float[N]; std::fill_n(direction,N,notanumber);
    Float* __restrict__ zirconium = new (std::nothrow) Float[N]; std::fill_n(zirconium,N,notanumber);
    Float* __restrict__ auxiliary = new (std::nothrow) Float[N]; std::fill_n(auxiliary,N,notanumber);
    Float* __restrict__ mittlerer = new (std::nothrow) Float[N]; std::fill_n(mittlerer,N,notanumber);
    assert( direction );
    assert( zirconium );
    assert( mittlerer );
    assert( auxiliary );
    
    for( int i = 0; i < N; i++ )
        direction[i] = zirconium[i] = auxiliary[i] = mittlerer[i] = 0.;
    
    Float z_r = notanumber;

    int K = 0;
    
    if( do_print_begin ) 
        LOGPRINTF( "(%d/%d)     BEGIN: Conjugate Gradient CSR (SSOR)\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csr_restart_on_full_dimension and K % N == 0 );
        
        bool preconresidual_seems_small = ( K != 0 ) and absolute(z_r) < tolerance*tolerance;

        if( restart_condition or ( preconresidual_seems_small and csr_restart_before_finish ) ) UNLIKELY {
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif
            for( int c = 0; c < N; c++ ) {
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
                residual[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
            }
            
            // inv( L^t + D/omega ) * D * inv( L + D/omega )
            
            // NOTE: Don't parallelize
            for( int c = 0; c < N; c++ ) {
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                Float aux = residual[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    if( csrcolumns[d] < c )
                        aux -= csrvalues[ d ] * mittlerer[ csrcolumns[d] ];
                
                mittlerer[c] = aux * omega /  diagonal[c];
                
            }
        
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif
            for( int c = 0; c < N; c++ ) {
            
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                mittlerer[c] *= diagonal[c];
                
                mittlerer[c] *= ( 2. - omega ) / omega;
            }
            
            // NOTE: Don't parallelize
            for( int c = N-1; c >= 0; c-- ) {
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
                Float aux = mittlerer[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    if( csrcolumns[d] > c )
                        aux -= csrvalues[ d ] * zirconium[ csrcolumns[d] ];
                
                zirconium[c] = aux * omega /  diagonal[c];
            
            }
            
            z_r = 0.;
        
            #if defined(_OPENMP)
            #pragma omp parallel for reduction( + : z_r )
            #endif
            for( int c = 0; c < N; c++ ) {
                            
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
                direction[c] = zirconium[c];
                       
                z_r += zirconium[c] * residual[c];
            
            }

            if( do_print_restart ) UNLIKELY 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );

            
        }
        
        /* printing information */

        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );

        if( print_condition and do_print_interim ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );
        }

        /* Check whether residual is small */
                
        bool preconresidual_is_small = absolute(z_r) < tolerance*tolerance;
        
        if( preconresidual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        // NOTE The calculation of d_r is reduced to r_r, which is already known.
        
        Float d_Ad = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction( + : d_Ad ) 
        #endif
        for( int c = 0; c < N; c++ )
        {
            auxiliary[c] = 0.;
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                auxiliary[c] += csrvalues[ d ] * direction[ csrcolumns[d] ];
                    
            d_Ad += direction[c] * auxiliary[c];
        }
        
        
        bool numerator_is_unreasonable   = not std::isfinite(z_r) or z_r < 0.;
        bool denominator_is_unreasonable = not std::isfinite(d_Ad) or d_Ad < 0.;
        bool denominator_is_small        = std::sqrt(absolute(d_Ad)) < machine_epsilon;
        
        if( numerator_is_unreasonable or denominator_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient energy is unreasonable with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        if( denominator_is_small ) UNLIKELY {
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient energy is small with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        
        const Float alpha = z_r / d_Ad;
    
        // NOTE: Don't parallelize 
        for( int c = 0; c < N; c++ )
        {
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            x[c] += alpha * direction[c];
            
            residual[c] -= alpha * auxiliary[c];
            
            Float aux = residual[c];
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                if( csrcolumns[d] < c )
                    aux -= csrvalues[ d ] * mittlerer[ csrcolumns[d] ];
            
            mittlerer[c] = aux * omega / diagonal[c];
        
        }
        
        Float z_r_new = 0.;
        
//         for( int c = 0; c < N; c++ ) {
//         
//         }
        
        // NOTE: Don't parallelize
        for( int c = N-1; c >= 0; c-- ) {
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            mittlerer[c] *= diagonal[c];
            
            mittlerer[c] *= ( 2. - omega ) / omega;
        
            Float aux = mittlerer[c];
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                if( csrcolumns[d] > c )
                    aux -= csrvalues[ d ] * zirconium[ csrcolumns[d] ];
            
            zirconium[c] = aux * omega / diagonal[c];
            
            z_r_new += zirconium[c] * residual[c];
            
        }
        
        const Float beta = z_r_new / z_r;
        
        z_r = z_r_new;
        
        #if defined(_OPENMP)
        #pragma omp parallel for 
        #endif
        for( int c = 0; c < N; c++ ) {
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            direction[c] = zirconium[c] + beta * direction[c];
            
        }
        
        K++;
        
    }
    
    if( do_print_finish ) UNLIKELY 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(z_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );

    
    delete[] ( direction ); 
    delete[] ( zirconium ); 
    delete[] ( mittlerer ); 
    delete[] ( auxiliary );

    return K;

}





int ConjugateGradientSolverCSR_SSOR_Eisenstat( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float precision,
    int print_modulo,
    const Float* __restrict__ diagonal,
    Float omega
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    assert( precision > 0 );
    // assert( print_modulo >= -1 );
    assert( diagonal );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = precision * sqrt( tolerance );
    tolerance = maximum( tolerance, desired_precision );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    const bool do_print_breakdown = print_modulo >= -2;
    const bool do_print_warning   = print_modulo >= -2;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    Float* __restrict__ direction = new (std::nothrow) Float[N]; std::fill_n(direction,N,notanumber);
    Float* __restrict__ zirconium = new (std::nothrow) Float[N]; std::fill_n(zirconium,N,notanumber);
    Float* __restrict__ technical = new (std::nothrow) Float[N]; std::fill_n(technical,N,notanumber);
    Float* __restrict__ auxiliary = new (std::nothrow) Float[N]; std::fill_n(auxiliary,N,notanumber);
    assert( direction );
    assert( zirconium );
    assert( technical );
    assert( auxiliary );
    
//     for( int i = 0; i < N; i++ )
//         direction[i] = zirconium[i] = technical[i] = auxiliary[i] = 0.; 
    
    Float z_r = notanumber;

    int K = 0;
    
    if( do_print_begin ) 
        LOGPRINTF( "(%d/%d)     BEGIN: Conjugate Gradient CSR (SSOR-Eisenstat)\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csr_restart_on_full_dimension and K % N == 0 );
        
        bool preconresidual_seems_small = ( K != 0 ) and absolute(z_r) < tolerance*tolerance;

        if( restart_condition or ( preconresidual_seems_small and csr_restart_before_finish ) ) UNLIKELY {
            

            // compute the residual, save it in 'auxiliary'
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif
            for( int c = 0; c < N; c++ ) {
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
                Float aux = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    if( diagonal[csrcolumns[d]] != 0. )
                        aux -= csrvalues[ d ] * x[ csrcolumns[d] ];

                auxiliary[c] = aux;
                
            }

            // Lower-left solve, save it in 'residual'
            // NOTE: Don't parallelize
            for( int c = 0; c < N; c++ ) {
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                    
                Float aux = auxiliary[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    if( csrcolumns[d] < c && diagonal[csrcolumns[d]] != 0. )
                        aux -= csrvalues[ d ] * residual[ csrcolumns[d] ];
                
                residual[c] = aux * omega / diagonal[c];
                
            }
        
            z_r = 0.;
        
            // compute z_r, z, and d from the vector r
            #if defined(_OPENMP)
            #pragma omp parallel for reduction( + : z_r )
            #endif
            for( int c = 0; c < N; c++ ) {
            
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                zirconium[c] = residual[c] * diagonal[c] * ( 2. - omega ) / omega;
                
                direction[c] = zirconium[c];
            
                z_r += zirconium[c] * residual[c];
            }
            
            if( do_print_restart ) UNLIKELY 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );

            
        }
        
        /* printing information */

        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );

        if( print_condition and do_print_interim ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );
        }

        /* Check whether residual is small */
                
        bool preconresidual_is_small = absolute(z_r) < tolerance*tolerance;
        
        if( preconresidual_is_small ) UNLIKELY 
            break;



        /* now the main work of the entire algorithm */
        
        // NOTE The calculation of d_r is reduced to z_r, which is already known.
        
        Float d_Ad = 0.;

        // NOTE: Don't parallelize
        // Upper right solve, save in 'technical'
        for( int c = N-1; c >= 0; c-- ) {
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            Float aux = direction[c];
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                if( csrcolumns[d] > c && diagonal[csrcolumns[d]] != 0. )
                    aux -= csrvalues[ d ] * technical[ csrcolumns[d] ];
            
            technical[c] = aux * omega / diagonal[c];
            
        }

        // NOTE: Don't parallelize (UPDATE d_Ad )
        // Lower left solve, save in 'auxiliary'
        for( int c = 0; c < N; c++ )
        {
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            Float aux = direction[c] - ( 2. / omega - 1. ) * diagonal[c] * technical[c];
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                if( csrcolumns[d] < c && diagonal[csrcolumns[d]] != 0. )
                    aux -= csrvalues[ d ] * auxiliary[ csrcolumns[d] ];
            
            auxiliary[c] = aux * omega / diagonal[c];
            
        }
        
        
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction( + : d_Ad )
        #endif
        for( int c = 0; c < N; c++ )
        {
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            auxiliary[c] += technical[c];

            d_Ad += auxiliary[c] * direction[c];

        } // NOTE: with extra memory or extra additions, we can get rid of the extra loop 
        
        
        
        bool numerator_is_unreasonable   = not std::isfinite(z_r) or z_r < 0.;
        bool denominator_is_unreasonable = not std::isfinite(d_Ad) or d_Ad < 0.;
        bool denominator_is_small        = std::sqrt(absolute(d_Ad)) < machine_epsilon;
        
        if( numerator_is_unreasonable or denominator_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient energy is unreasonable with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        if( denominator_is_small and do_print_warning ) UNLIKELY {
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient energy is small with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        
        const Float alpha = z_r / d_Ad;
    
        Float z_r_new = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction( + : z_r_new )
        #endif
        for( int c = 0; c < N; c++ )
        {
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            x[c] += alpha * technical[c];
            
            residual[c] -= alpha * auxiliary[c];

            zirconium[c] = residual[c] * ( diagonal[c] * ( 2. - omega ) / omega );

            z_r_new += zirconium[c] * residual[c];
            
        }
        
        const Float beta = z_r_new / z_r;
        
        z_r = z_r_new;
        
        #if defined(_OPENMP)
        #pragma omp parallel for 
        #endif
        for( int c = 0; c < N; c++ ) {
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            direction[c] = zirconium[c] + beta * direction[c];
            
        }
        
        K++;
        
    }
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(z_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );

    
    delete[] ( direction ); 
    delete[] ( zirconium ); 
    delete[] ( auxiliary );
    delete[] ( technical );

    return K;

}






int ConjugateGradientSolverCSR_Rainbow( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float precision,
    int print_modulo,
    const Float* __restrict__ diagonal,
    Float omega,
    int num_colors, const int* __restrict__ F, const int* __restrict__ B, const int* __restrict__ R
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    assert( precision > 0 );
    // assert( print_modulo >= -1 );
    assert( diagonal );

    assert( num_colors >= 0 );
    assert( F );
    assert( B );
    assert( R );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = precision * sqrt( tolerance );
    tolerance = maximum( tolerance, desired_precision );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    const bool do_print_breakdown = print_modulo >= -2;
    const bool do_print_warning   = print_modulo >= -2;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    Float* __restrict__ direction = new (std::nothrow) Float[N]; std::fill_n(direction,N,notanumber);
    Float* __restrict__ zirconium = new (std::nothrow) Float[N]; std::fill_n(zirconium,N,notanumber);
    Float* __restrict__ auxiliary = new (std::nothrow) Float[N]; std::fill_n(auxiliary,N,notanumber);
    Float* __restrict__ mittlerer = new (std::nothrow) Float[N]; std::fill_n(mittlerer,N,notanumber);
    assert( direction );
    assert( zirconium );
    assert( mittlerer );
    assert( auxiliary );
    
    for( int i = 0; i < N; i++ )
        direction[i] = zirconium[i] = auxiliary[i] = mittlerer[i] = 0.;
    
    Float z_r = notanumber;

    int K = 0;
    
    if( do_print_begin ) 
        LOGPRINTF( "(%d/%d)     BEGIN: Conjugate Gradient CSR (SSOR-Rainbow)\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csr_restart_on_full_dimension and K % N == 0 );
        
        bool preconresidual_seems_small = ( K != 0 ) and absolute(z_r) < tolerance*tolerance;

        if( restart_condition or ( preconresidual_seems_small and csr_restart_before_finish ) ) UNLIKELY {
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif
            for( int c = 0; c < N; c++ ) {
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
                residual[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
            }
            
            // (2-omega)/omega * inv( L^t + D/omega ) * D * inv( L + D/omega )
            
            // NOTE: Rainbow-ing of CSR matrix used here 
            for( int color = 0; color < num_colors; color++ ) {
            
                #if defined(_OPENMP)
                #pragma omp parallel for
                #endif
                for( int i = B[color]; i < B[color+1]; i++ ) {

                    int c = R[i];
                    
                    // assert( F[c] == color );
                    
                    if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 

                    Assert( F[c] == color, N, c, F[c], color );
                    
                    Float aux = residual[c];
                    
                    for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                        if( F[csrcolumns[d]] < F[c] )
                            aux -= csrvalues[ d ] * mittlerer[ csrcolumns[d] ];
                    
                    mittlerer[c] = aux * omega /  diagonal[c];
                    
                }

            }
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif
            for( int c = 0; c < N; c++ ) {
            
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                mittlerer[c] *= diagonal[c];
                
                mittlerer[c] *= ( 2. - omega ) / omega;
            }
            
            // NOTE: Rainbow-ing of CSR matrix used here 
            for( int color = num_colors - 1; color >= 0; color-- ) {
            
                #if defined(_OPENMP)
                #pragma omp parallel for
                #endif
                for( int i = B[color]; i < B[color+1]; i++ ) {

                    int c = R[i];
                    
                    if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                    // assert( F[c] == color );
                    
                    Float aux = mittlerer[c];
                    
                    for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                        if( F[csrcolumns[d]] > F[c] )
                            aux -= csrvalues[ d ] * zirconium[ csrcolumns[d] ];
                    
                    zirconium[c] = aux * omega /  diagonal[c];
                
                }
            
            }
            
            z_r = 0.;
        
            #if defined(_OPENMP)
            #pragma omp parallel for reduction( + : z_r )
            #endif
            for( int c = 0; c < N; c++ ) {
                            
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
                direction[c] = zirconium[c];
                       
                z_r += zirconium[c] * residual[c];
            
            }

            if( do_print_restart ) UNLIKELY 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );

            
        }
        
        /* printing information */

        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );

        if( print_condition and do_print_interim ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );
        }
        
        /* Check whether residual is small */
                
        bool preconresidual_is_small = absolute(z_r) < tolerance*tolerance;
        
        if( preconresidual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        // NOTE The calculation of d_r is reduced to r_r, which is already known.
        
        Float d_Ad = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction( + : d_Ad ) 
        #endif
        for( int c = 0; c < N; c++ )
        {
            auxiliary[c] = 0.;
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                auxiliary[c] += csrvalues[ d ] * direction[ csrcolumns[d] ];
                    
            d_Ad += direction[c] * auxiliary[c];
        }
        
        
        bool numerator_is_unreasonable   = not std::isfinite(z_r) or z_r < 0.;
        bool denominator_is_unreasonable = not std::isfinite(d_Ad) or d_Ad < 0.;
        bool denominator_is_small        = std::sqrt(absolute(d_Ad)) < machine_epsilon;
        
        if( numerator_is_unreasonable or denominator_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient energy is unreasonable with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        if( denominator_is_small ) UNLIKELY {
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient energy is small with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        
        const Float alpha = z_r / d_Ad;
    
        // NOTE: Rainbow-ing of CSR matrix used here 
        for( int color = 0; color < num_colors; color++ ) {
        
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif
            for( int i = B[color]; i < B[color+1]; i++ ) {

                int c = R[i];
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                // assert( F[c] == color );
                    
                x[c] += alpha * direction[c];
                
                residual[c] -= alpha * auxiliary[c];
                
                Float aux = residual[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    if( F[csrcolumns[d]] < F[c] )
                        aux -= csrvalues[ d ] * mittlerer[ csrcolumns[d] ];
                
                mittlerer[c] = aux * omega / diagonal[c];
            
            }

        }
        
        Float z_r_new = 0.;
        
//         for( int c = 0; c < N; c++ ) {
//         
//         }
        
        // NOTE: Rainbow-ing of CSR matrix used here 
        for( int color = num_colors - 1; color >= 0; color-- ) {
        
            Float z_r_new_local = 0.; 

            #if defined(_OPENMP)
            #pragma omp parallel for reduction( + : z_r_new_local ) 
            #endif
            for( int i = B[color]; i < B[color+1]; i++ ) {

                int c = R[i];
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                // assert( F[c] == color );
                    
                mittlerer[c] *= diagonal[c];
                
                mittlerer[c] *= ( 2. - omega ) / omega;
            
                Float aux = mittlerer[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    if( F[csrcolumns[d]] > F[c] )
                        aux -= csrvalues[ d ] * zirconium[ csrcolumns[d] ];
                
                zirconium[c] = aux * omega / diagonal[c];
                
                z_r_new_local += zirconium[c] * residual[c];
                
            }

            z_r_new += z_r_new_local;

        }
        
        const Float beta = z_r_new / z_r;
        
        z_r = z_r_new;
        
        #if defined(_OPENMP)
        #pragma omp parallel for 
        #endif
        for( int c = 0; c < N; c++ ) {
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            direction[c] = zirconium[c] + beta * direction[c];
            
        }
        
        K++;
        
    }
    
    if( do_print_finish ) UNLIKELY 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(z_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );

    
    delete[] ( direction ); 
    delete[] ( zirconium ); 
    delete[] ( mittlerer ); 
    delete[] ( auxiliary );

    return K;

}






int ConjugateGradientSolverCSR_Eisenstat_Rainbow( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float precision,
    int print_modulo,
    const Float* __restrict__ diagonal,
    Float omega,
    int num_colors, const int* __restrict__ F, const int* __restrict__ B, const int* __restrict__ R
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    assert( precision > 0 );
    // assert( print_modulo >= -1 );
    assert( diagonal );

    assert( num_colors >= 0 );
    assert( F );
    assert( B );
    assert( R );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = precision * sqrt( tolerance );
    tolerance = maximum( tolerance, desired_precision );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    const bool do_print_breakdown = print_modulo >= -2;
    const bool do_print_warning   = print_modulo >= -2;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    Float* __restrict__ direction = new (std::nothrow) Float[N]; std::fill_n(direction,N,notanumber);
    Float* __restrict__ zirconium = new (std::nothrow) Float[N]; std::fill_n(zirconium,N,notanumber);
    Float* __restrict__ technical = new (std::nothrow) Float[N]; std::fill_n(technical,N,notanumber);
    Float* __restrict__ auxiliary = new (std::nothrow) Float[N]; std::fill_n(auxiliary,N,notanumber);
    assert( direction );
    assert( zirconium );
    assert( technical );
    assert( auxiliary );
    
    Float z_r = notanumber;

    for( int c = 0; c < N; c++ ) direction[c] = zirconium[c] = technical[c] = auxiliary[c] = notanumber;

    int K = 0;
    
    if( do_print_begin ) 
        LOGPRINTF( "(%d/%d)     BEGIN: Conjugate Gradient CSR (SSOR-Rainbow-Eisenstat)\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csr_restart_on_full_dimension and K % N == 0 );
        
        bool preconresidual_seems_small = ( K != 0 ) and absolute(z_r) < tolerance*tolerance;

        if( restart_condition or ( preconresidual_seems_small and csr_restart_before_finish ) ) UNLIKELY {
            

            // compute the residual, save it in 'auxiliary'
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif
            for( int c = 0; c < N; c++ ) {
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
                Float aux = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    if( diagonal[csrcolumns[d]] != 0. )
                        aux -= csrvalues[ d ] * x[ csrcolumns[d] ];

                auxiliary[c] = aux;
                
            }

            // Lower-left solve, from 'auxiliary' into 'residual'
            // NOTE: Rainbow-ing of CSR matrix used here 
            for( int color = 0; color < num_colors; color++ ) {
            
                #if defined(_OPENMP)
                #pragma omp parallel for
                #endif
                for( int i = B[color]; i < B[color+1]; i++ ) {

                    const int c = R[i]; // assert( F[c] == color );
                    
                    if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                    
                    Float aux = auxiliary[c];
                    
                    for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                        if( F[csrcolumns[d]] < F[c] && diagonal[csrcolumns[d]] != 0. )
                            aux -= csrvalues[ d ] * residual[ csrcolumns[d] ];
                    
                    residual[c] = aux * omega / diagonal[c];

                    assert( std::isfinite(residual[c]) );
                
                }

            }

            z_r = 0.;
        
            // compute z_r, z, and d from the vector r
            #if defined(_OPENMP)
            #pragma omp parallel for reduction( + : z_r )
            #endif
            for( int c = 0; c < N; c++ ) {
            
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                zirconium[c] = residual[c] * diagonal[c] * ( 2. - omega ) / omega;
                
                direction[c] = zirconium[c];
            
                z_r += zirconium[c] * residual[c];
            }
            
            if( do_print_restart ) 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );

            
        }
        
        /* printing information */

        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );

        if( print_condition and do_print_interim ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );
        }

        /* Check whether residual is small */
                
        bool preconresidual_is_small = absolute(z_r) < tolerance*tolerance;
        
        if( preconresidual_is_small )
            break;



        /* now the main work of the entire algorithm */
        
        // NOTE The calculation of d_r is reduced to z_r, which is already known.
        
        Float d_Ad = 0.;

        // Upper right solve, from 'direction' into 'technical'
        // NOTE: Rainbow-ing of CSR matrix used here 
        for( int color = num_colors - 1; color >= 0; color-- ) {
        
            // Float z_r_new_local = 0.; // TODO(martinlicht): understand whether that variable is needed or whether it is a relic

            #if defined(_OPENMP)
            // #pragma omp parallel for reduction( + : z_r_new_local ) 
            #endif
            for( int i = B[color]; i < B[color+1]; i++ ) {

                int c = R[i]; // assert( F[c] == color );
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                Float aux = direction[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    if( F[csrcolumns[d]] > F[c] )
                        aux -= csrvalues[ d ] * technical[ csrcolumns[d] ];
                
                technical[c] = aux * omega / diagonal[c];

                assert( std::isfinite(technical[c]) );
                
            }

        }
        
        
        // Lower left solve, into 'auxiliary'
        // NOTE: Rainbow-ing of CSR matrix used here 
        for( int color = 0; color < num_colors; color++ ) {
        
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif
            for( int i = B[color]; i < B[color+1]; i++ ) {

                int c = R[i]; // assert( F[c] == color );
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                Float aux = direction[c] - ( 2. / omega - 1. ) * diagonal[c] * technical[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    if( F[csrcolumns[d]] < F[c] )
                        aux -= csrvalues[ d ] * auxiliary[ csrcolumns[d] ];
                
                auxiliary[c] = aux * omega / diagonal[c];
            
            }

        }
        
        
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction( + : d_Ad )
        #endif
        for( int c = 0; c < N; c++ )
        {
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 

            auxiliary[c] += technical[c];

            d_Ad += auxiliary[c] * direction[c];

        } // NOTE: with extra memory or extra additions, we can get rid of the extra loop 
        
        
        
        bool numerator_is_unreasonable   = not std::isfinite(z_r) or z_r < 0.;
        bool denominator_is_unreasonable = not std::isfinite(d_Ad) or d_Ad < 0.;
        bool denominator_is_small        = std::sqrt(absolute(d_Ad)) < machine_epsilon;
        
        if( numerator_is_unreasonable or denominator_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient energy is unreasonable with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        if( denominator_is_small ) UNLIKELY {
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient energy is small with %.9le\n", K, N, (double)(safedouble)d_Ad );
            break;
        }
        
        
        const Float alpha = z_r / d_Ad;
    
        Float z_r_new = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction( + : z_r_new )
        #endif
        for( int c = 0; c < N; c++ )
        {
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            x[c] += alpha * technical[c];
            
            residual[c] -= alpha * auxiliary[c];

            zirconium[c] = residual[c] * ( diagonal[c] * ( 2. - omega ) / omega );

            z_r_new += zirconium[c] * residual[c];
            
        }
        
        const Float beta = z_r_new / z_r;
        
        z_r = z_r_new;
        
        #if defined(_OPENMP)
        #pragma omp parallel for 
        #endif
        for( int c = 0; c < N; c++ ) {
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            direction[c] = zirconium[c] + beta * direction[c];
            
        }
        
        K++;
        
    }
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(z_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble) std::sqrt(z_r), (double)(safedouble) tolerance );

    
    delete[] ( direction ); 
    delete[] ( zirconium ); 
    delete[] ( auxiliary );
    delete[] ( technical );

    return K;

}









































int ConjugateResidualSolverCSR( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ res,
    const Float precision,
    int print_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( res );
    assert( precision > 0 );
    // assert( print_modulo >= -1 );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = precision * sqrt( tolerance );
    tolerance = maximum( tolerance, desired_precision );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    const bool do_print_breakdown = print_modulo >= -2;
    const bool do_print_warning   = print_modulo >= -2;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    Float* __restrict__  dir = new (std::nothrow) Float[N]; std::fill_n( dir,N,notanumber);
    Float* __restrict__ Adir = new (std::nothrow) Float[N]; std::fill_n(Adir,N,notanumber);
    Float* __restrict__ Ares = new (std::nothrow) Float[N]; std::fill_n(Ares,N,notanumber);
    assert(  dir );
    assert( Adir );
    assert( Ares );
    
    Float* __restrict__  vil = new (std::nothrow) Float[N];
    assert( vil );
    
    
    Float Ad_r  = notanumber;
    Float Ad_Ad = notanumber;
    
    int K = 0;
    
    if( do_print_begin ) 
        LOGPRINTF( "(%d/%d)     BEGIN: Conjugate Residual CSR\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csr_restart_on_full_dimension and K % N == 0 );
        
        bool residualenergy_seems_small = ( K != 0 ) and absolute(Ad_r) < tolerance*tolerance;

        if( restart_condition or ( residualenergy_seems_small and csr_restart_before_finish ) ) UNLIKELY {
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                res[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    res[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                dir[c] = res[c]; // this line seems to slow down performance ....
                
            }
            
            Ad_r  = 0.;
            Ad_Ad = 0.;
            
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(+:Ad_r,Ad_Ad)
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                Adir[c] = 0.;
                
                for( int e = csrrows[c]; e < csrrows[c+1]; e++ )
                    Adir[c] += csrvalues[ e ] * dir[ csrcolumns[e] ];
                
                Ares[c] = Adir[c];
                
                Ad_r  += Adir[c] *  res[c];
                Ad_Ad += Adir[c] * Adir[c];
                
            }
            
            if( do_print_restart ) 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Ad_r), (double)(safedouble) tolerance );
            
        }
        
        /* printing information */

        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );

        if( print_condition and do_print_interim ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Ad_r), (double)(safedouble) tolerance );
        }
        
        /* Check whether res is small */
                
        bool residualenergy_is_small = absolute(Ad_r) < tolerance*tolerance;
        bool residualenergy_is_unreasonable = not std::isfinite(Ad_r) or Ad_r < 0.;
        
        if( residualenergy_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Residual energy is unreasonable with %.9le\n", K, N, (double)(safedouble)Ad_r );
            break;
        }

        if( residualenergy_is_small )
            break;
        
        
        
        bool denominator_is_unreasonable = not std::isfinite(Ad_Ad) or Ad_Ad < 0.;
        bool denominator_is_small        = std::sqrt(absolute(Ad_Ad)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient double energy is unreasonable with %.9le\n", K, N, (double)(safedouble)Ad_Ad );
            break;
        }
        
        if( denominator_is_small ) UNLIKELY {
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Ad_r), (double)(safedouble) tolerance );
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient double energy is small with %.9le\n", K, N, (double)(safedouble)Ad_Ad );
            break;
        }
        
        
        

        /* now the main work of the entire algorithm */
        
        Float alpha = Ad_r / Ad_Ad;
        
        Float new_Ar_r = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:new_Ar_r)
        #endif 
        for( int c = 0; c < N; c++ )
        {
            vil[c] = 0.;
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                vil[c] += csrvalues[ d ] * Adir[ csrcolumns[d] ];
                    
            x[c] = x[c] + alpha * dir[c];
            
            res[c] = res[c] - alpha * Adir[c];
            
            Ares[c] = Ares[c] - alpha * vil[c];
            
            new_Ar_r += Ares[c] * res[c];
            
        }
        
        Float beta = new_Ar_r / Ad_r;
        
        Ad_Ad = 0.;
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:Ad_Ad)
        #endif
        for( int c = 0; c < N; c++ )
        {
            
            dir[c]  =  res[c] + beta *  dir[c];
            
            Adir[c] = Ares[c] + beta * Adir[c];
            
            Ad_Ad += Adir[c] * Adir[c];
        }
        
        Ad_r = new_Ar_r;
                
        K++;
        
    }
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(Ad_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble) std::sqrt(Ad_r), (double)(safedouble) tolerance );

    
    delete[] (  dir );
    delete[] ( Adir );
    delete[] ( Ares );

    delete[] ( vil );

    return K;

}




























int ConjugateResidualSolverCSR_textbook( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ res,
    const Float precision,
    int print_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( res );
    assert( precision > 0 );
    // assert( print_modulo >= -1 );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = precision * sqrt( tolerance );
    tolerance = maximum( tolerance, desired_precision );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    const bool do_print_breakdown = print_modulo >= -2;
    const bool do_print_warning   = print_modulo >= -2;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    Float* __restrict__  dir = new (std::nothrow) Float[N]; std::fill_n( dir,N,notanumber);
    Float* __restrict__ Adir = new (std::nothrow) Float[N]; std::fill_n(Adir,N,notanumber);
    Float* __restrict__ Ares = new (std::nothrow) Float[N]; std::fill_n(Ares,N,notanumber);
    assert(  dir );
    assert( Adir );
    assert( Ares );
    
    Float* __restrict__  vil = new (std::nothrow) Float[N];
    assert( vil );
    
    
    Float Ar_r  = notanumber;
    Float Ad_Ad = notanumber;

    int K = 0;
    
    if( do_print_begin ) 
        LOGPRINTF( "(%d/%d)     BEGIN: Conjugate Residual CSR (textbook)\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csr_restart_on_full_dimension and K % N == 0 );
        
        bool residualenergy_seems_small = ( K != 0 ) and absolute(Ar_r) < tolerance*tolerance;

        if( restart_condition or ( residualenergy_seems_small and csr_restart_before_finish ) ) UNLIKELY {
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                res[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    res[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                dir[c] = res[c]; // this line seems to slow down performance ....
                
            }
            
            Ar_r  = 0.;
            Ad_Ad = 0.;
            
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(+:Ar_r,Ad_Ad)
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                Adir[c] = 0.;
                
                for( int e = csrrows[c]; e < csrrows[c+1]; e++ )
                    Adir[c] += csrvalues[ e ] * dir[ csrcolumns[e] ];
                
                Ares[c] = Adir[c];
                
                Ar_r  += Ares[c] *  res[c];
                Ad_Ad += Adir[c] * Adir[c];
                
            }

            if( do_print_restart ) 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Ar_r), (double)(safedouble) tolerance );
            
        }
        
        /* printing information */

        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );

        if( print_condition and do_print_interim ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Ar_r), (double)(safedouble) tolerance );
        }

        /* Check whether res is small */
                
        bool residualenergy_is_unreasonable = not std::isfinite(Ar_r) or Ar_r < 0.;
        bool residualenergy_is_small = absolute(Ar_r) < tolerance*tolerance;
        
        if( residualenergy_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Residual energy is unreasonable with %.9le\n", K, N, (double)(safedouble)Ar_r );
            break;
        }
  
        if( residualenergy_is_small )
            break;


        bool denominator_is_unreasonable = not std::isfinite(Ad_Ad) or Ad_Ad < 0.;
        bool denominator_is_small        = std::sqrt(absolute(Ad_Ad)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) UNLIKELY {
            if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient double energy is unreasonable with %.9le\n", K, N, (double)(safedouble)Ad_Ad );
            break;
        }
        
        if( denominator_is_small ) UNLIKELY {
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) std::sqrt(Ar_r), (double)(safedouble) tolerance );
            if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient double energy is small with %.9le\n", K, N, (double)(safedouble)Ad_Ad );
            break;
        }
        
        
        /* now the main work of the entire algorithm */
        
        Float alpha = Ar_r / Ad_Ad;
        
        Float new_Ar_r = 0.;
        
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:new_Ar_r)
        #endif 
        for( int c = 0; c < N; c++ )
        {
            vil[c] = 0.;
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                vil[c] += csrvalues[ d ] * Adir[ csrcolumns[d] ];
                    
            x[c] = x[c] + alpha * dir[c];
            
            res[c] = res[c] - alpha * Adir[c];
            
            Ares[c] = Ares[c] - alpha * vil[c];
            
            new_Ar_r += Ares[c] * res[c];
            
        }
        
        Float beta = new_Ar_r / Ar_r;
        
        Ad_Ad = 0.;
        #if defined(_OPENMP)
        #pragma omp parallel for reduction(+:Ad_Ad)
        #endif 
        for( int c = 0; c < N; c++ )
        {
            
            dir[c]  =  res[c] + beta *  dir[c];
            
            Adir[c] = Ares[c] + beta * Adir[c];
            
            Ad_Ad += Adir[c] * Adir[c];
        }
        
        Ar_r = new_Ar_r;
                
        K++;
        
    }
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(Ar_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble) std::sqrt(Ar_r), (double)(safedouble) tolerance );

    
    delete[] (  dir );
    delete[] ( Adir );
    delete[] ( Ares );

    delete[] ( vil );

    return K;

}




































                      




int MINRESCSR( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ res,
    const Float precision,
    int print_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( res );
    assert( precision > 0 );
    // assert( print_modulo >= -1 );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = precision * sqrt( tolerance );
    tolerance = maximum( tolerance, desired_precision );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    // const bool do_print_breakdown = print_modulo >= -2;
    // const bool do_print_warning   = print_modulo >= -2;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    Float* __restrict__ v0 = new (std::nothrow) Float[N]; std::fill_n(v0,N,notanumber);
    Float* __restrict__ v1 = new (std::nothrow) Float[N]; std::fill_n(v1,N,notanumber);
    Float* __restrict__ w0 = new (std::nothrow) Float[N]; std::fill_n(w0,N,notanumber);
    Float* __restrict__ w1 = new (std::nothrow) Float[N]; std::fill_n(w1,N,notanumber);
    
    assert( v0 );
    assert( v1 );
    assert( w0 );
    assert( w1 );
    
    Float* __restrict__ vn = new (std::nothrow) Float[N]; std::fill_n(vn,N,notanumber);
    Float* __restrict__ wn = new (std::nothrow) Float[N]; std::fill_n(wn,N,notanumber);
    Float* __restrict__  p = new (std::nothrow) Float[N]; std::fill_n( p,N,notanumber);
    
    assert( vn );
    assert( wn );
    assert(  p );
    
    Float gamma = notanumber;
    Float eta   = notanumber;
    
    Float s0 = notanumber;
    Float s1 = notanumber;
    Float c0 = notanumber;
    Float c1 = notanumber;

    int K = 0;

    if( do_print_begin ) 
        LOGPRINTF( "(%d/%d)     BEGIN: MINRES CSR\n", K, N );

    while( K < N ){
        
        
        bool restart_condition = (K == 0) or ( csr_restart_on_full_dimension and K % N == 0 );
        
        bool residual_seems_small = ( K != 0 ) and ( absolute(eta) < tolerance );
        
        if( restart_condition or ( residual_seems_small and csr_restart_before_finish ) ) {
            
            Float gamma_sq = 0.;
            
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(+:gamma_sq)
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                v0[c] = 0.;
                w0[c] = 0.;
                v1[c] = b[c];
                w1[c] = 0.;
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    v1[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                gamma_sq += v1[c] * v1[c];
                
            }
            
            assert( gamma_sq > 0. );
            
            gamma = std::sqrt(gamma_sq);
            
            #if defined(_OPENMP)
            #pragma omp parallel for 
            #endif 
            for( int c = 0; c < N; c++ ) 
                v1[c] /= gamma;
            
            s0 = s1 = 0.;
            c0 = c1 = 1.;
            
            eta = gamma;

            if( do_print_restart ) 
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble) absolute(eta), (double)(safedouble) tolerance );

            
        }
        
        bool residual_is_small = ( absolute(eta) < tolerance);
        
        if( residual_is_small )
            break;

//         if( K % 100 == 0 )
//             LOG << K << space << gamma << space << eta << space << gamma/eta << nl;
            
        {
            
            
            Float delta = 0.;
            
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(+:delta)
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                p[c] = 0.;
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    p[c] += csrvalues[ d ] * v1[ csrcolumns[d] ];
                
                delta += p[c] * v1[c];
                
            }
            
            assert( delta > 0. );
            
            Float gamma_n_sq = 0.;
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(+:gamma_n_sq)
            #endif 
            for( int c = 0; c < N; c++ ) {
                
                vn[c] = p[c] - delta * v1[c] - gamma * v0[c];
                
                gamma_n_sq += vn[c] * vn[c];
                
            }
            
            assert( gamma_n_sq > 0. );
 
            Float gamma_n = std::sqrt(gamma_n_sq);
            
            #if defined(_OPENMP)
            #pragma omp parallel for 
            #endif 
            for( int c = 0; c < N; c++ )
                vn[c] /= gamma_n;
            
            Float alpha_0 = c1 * delta - c0 * s1 * gamma;
            Float alpha_1 = std::sqrt( alpha_0 * alpha_0 + gamma_n * gamma_n );
            Float alpha_2 = s1 * delta + c0 * c1 * gamma;
            Float alpha_3 = s0 * gamma;
 
            assert( alpha_1 > 0. );

            Float cn = alpha_0 / alpha_1;
            Float sn = gamma_n / alpha_1;
            
            #if defined(_OPENMP)
            #pragma omp parallel for
            #endif 
            for( int c = 0; c < N; c++ ) {
                wn[c] = ( v1[c] - alpha_2 * w1[c] - alpha_3 * w0[c] ) / alpha_1;
                x[c] = x[c] + cn * eta * wn[c];
            }
            
            eta = - sn * eta;
            
//             LOG << "\t" << alpha_0 << space << alpha_1 << space << alpha_2 << space << alpha_3 << nl;
//             LOG << "\t" << cn << space << sn << space << eta << nl;
            
            
            std::swap( v0, v1 );
            std::swap( w0, w1 );
            std::swap( v1, vn );
            std::swap( w1, wn );
//             v0 = v1;
//             w0 = w1;
//             v1 = vn;
//             w1 = wn;
            
            gamma = gamma_n;
            
            c0 = c1;
            s0 = s1;
            c1 = cn;
            s1 = sn;

        }

        
        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );

        if( print_condition and do_print_interim ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble)eta, (double)(safedouble) tolerance );
        }

        K++;
        
    }
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, eta < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble)eta, (double)(safedouble) tolerance );


    for( int c = 0; c < N; c++ ) {
        res[c] = b[c];
        for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
            res[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
    }

    
    delete[] ( v0 );
    delete[] ( v1 );
    delete[] ( w0 );
    delete[] ( w1 );
    
    delete[] ( vn );
    delete[] ( wn );
    delete[] (  p );

    return K;
    
}
























// The Convergence of Inexact Chebyshev and Richardson Iterative Methods for Solving Linear Systems


int ChebyshevIteration_DiagonalPreconditioner( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float precision,
    int print_modulo,
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
    assert( precision > 0 );
    assert( precon );
    
    Float tolerance = 0.;
    for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
    tolerance = precision * sqrt( tolerance );
    tolerance = maximum( tolerance, desired_precision );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    // const bool do_print_restart   = print_modulo >=  0;
    // const bool do_print_breakdown = print_modulo >= -2;
    // const bool do_print_warning   = print_modulo >= -2;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    Float* __restrict__ zaratite = new (std::nothrow) Float[N];
    assert( zaratite );
    
    Float* x_prev = new (std::nothrow) Float[N]; std::fill_n(x_prev,N,notanumber);
    Float* x_curr = new (std::nothrow) Float[N]; std::fill_n(x_curr,N,notanumber);
    Float* x_next = new (std::nothrow) Float[N]; std::fill_n(x_next,N,notanumber);
    assert( x_prev );
    assert( x_curr );
    assert( x_next );

    #if defined(_OPENMP)
    #pragma omp parallel for 
    #endif
    for( int c = 0; c < N; c++ )
        x_prev[c] = x_curr[c] = x_next[c] = x[c];
            
            

    LOG << lower << space << upper << nl;
    
    const Float alpha = 2. / ( upper + lower );
    const Float mu    = ( upper + lower ) / ( upper - lower );
    
    Float gamma_prev = notanumber;
    Float gamma_curr = notanumber;
    //Float gamma_next = notanumber;
    
    Float r_r = notanumber;
    
    int K = 0;
    
    if( do_print_begin ) 
        LOGPRINTF( "(%d/%d)     BEGIN: Chebyshev (Diag) CSR\n", K, N );

    while( K < N ){
        
        bool restart_condition = ( K == 0 ) or ( csr_restart_on_full_dimension and K % N == 0 );

        bool residual_seems_small = ( K != 0 ) and std::sqrt(r_r) < tolerance;

        if( restart_condition or ( residual_seems_small and csr_restart_before_finish ) ) UNLIKELY {
            
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
                
        bool residual_is_small = std::sqrt(r_r) < tolerance;
        
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
        
        gamma_prev = gamma_curr; 
        gamma_curr = gamma_next; 
        // std::swap( gamma_curr, gamma_prev );
        // std::swap( gamma_next, gamma_curr );
        std::swap( x_curr, x_prev );
        std::swap( x_next, x_curr );
        
        
        bool print_condition = ( print_modulo > 0 and K % print_modulo == 0 );

        if( print_condition and do_print_interim ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", K, N, (double)(safedouble)std::sqrt(r_r), (double)(safedouble) tolerance );
        }

        K++;
        
    }
    
    if( do_print_finish )
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", K, N, std::sqrt(r_r) < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble)std::sqrt(r_r), (double)(safedouble) tolerance );

    
    #if defined(_OPENMP)
    #pragma omp parallel for 
    #endif
    for( int c = 0; c < N; c++ )
        x[c] = x_next[c];
            
            
    delete[] ( x_prev );
    delete[] ( x_curr );
    delete[] ( x_next );

    delete[] ( zaratite );

    return K;

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
//     const Float precision = tolerance;
//     
//     /*************/
//     /* MAIN LOOP */
//     /*************/
//     
//     while( iter < maxiter ) {
//         
//         LOGPRINTF( "@@@@@@@@@@ Uzawa-CRM Iteration %d / %d: unorm is %f\n", iter, maxiter, vectornorm(pdim,u) );
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
//         LOGPRINTF( "@@@@@@@@@@ Residual norm: %f\n", residualnorm );
//         
//         if( residualnorm < tolerance ) {
//             LOGPRINTF( "@@@@@@@@@@ Threshold deceeded.\n" );
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
// //     LOGPRINTF( "@@@@@@@@@@ Residual sigma %f \n", vectornorm( M, ressigma ) );
// //     LOGPRINTF( "@@@@@@@@@@ Residual u     %f \n", vectornorm( pdim, resu  ) );
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
// //         LOGPRINTF( "@@@@@@@@@@ Other residual sigma %f \n", vectornorm( M, ressigma ) );
// //         LOGPRINTF( "@@@@@@@@@@ Other residual u     %f \n", vectornorm( pdim, resu  ) );
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
//         LOGPRINTF( "@@@@@@@@@@ One more Residual sigma %f \n", vectornorm( M, ressigma ) );
//         LOGPRINTF( "@@@@@@@@@@ One more Residual u     %f \n", vectornorm( pdim, resu  ) );
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
//     delete[] ( r       );
//     delete[] ( C_r     ); 
//     delete[] ( AiBt_r  );
//     delete[] ( BAiBt_r );
//     delete[] ( d       );
//     delete[] ( C_d     );
//     delete[] ( AiBt_d  );
//     delete[] ( BAiBt_d );
//     delete[] ( p       );
//     delete[] ( C_p     );
//     delete[] ( AiBt_p  );
//     delete[] ( BAiBt_p );
//     delete[] ( tempM );
//     delete[] ( tempN );
//     delete[] ( entriesBt );
//     delete[] ( csrrowsBt );
//     delete[] ( csrcolumnsBt );
//     
//     LOGPRINTF( "@@@@@@@@@@ Exit: %f %f\n", vectornorm( M, ressigma ), vectornorm( pdim, resu  ) );
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
//     const Float precision = tolerance;
//     
//     /*************/
//     /* MAIN LOOP */
//     /*************/
//     
//     while( iter < maxiter ) {
//         
//         LOGPRINTF( "@@@@@@@@@@ Uzawa-CRM Iteration %d / %d: unorm is %f\n", iter, maxiter, vectornorm(pdim,u) );
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
//         LOGPRINTF( "@@@@@@@@@@ Residual norm: %f\n", residualnorm );
//         
//         if( residualnorm < tolerance ) {
//             LOGPRINTF( "@@@@@@@@@@ Threshold deceeded.\n" );
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
// //     LOGPRINTF( "@@@@@@@@@@ Residual sigma %f \n", vectornorm( M, ressigma ) );
// //     LOGPRINTF( "@@@@@@@@@@ Residual u     %f \n", vectornorm( pdim, resu  ) );
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
// //         LOGPRINTF( "@@@@@@@@@@ Other residual sigma %f \n", vectornorm( M, ressigma ) );
// //         LOGPRINTF( "@@@@@@@@@@ Other residual u     %f \n", vectornorm( pdim, resu  ) );
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
//         LOGPRINTF( "@@@@@@@@@@ One more Residual sigma %f \n", vectornorm( M, ressigma ) );
//         LOGPRINTF( "@@@@@@@@@@ One more Residual u     %f \n", vectornorm( pdim, resu  ) );
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
//     delete[] ( r       );
//     delete[] ( C_r     ); 
//     delete[] ( AiBt_r  );
//     delete[] ( BAiBt_r );
//     delete[] ( d       );
//     delete[] ( C_d     );
//     delete[] ( AiBt_d  );
//     delete[] ( BAiBt_d );
//     delete[] ( p       );
//     delete[] ( C_p     );
//     delete[] ( AiBt_p  );
//     delete[] ( BAiBt_p );
//     delete[] ( tempM );
//     delete[] ( tempN );
//     delete[] ( entriesBt );
//     delete[] ( csrrowsBt );
//     delete[] ( csrcolumnsBt );
//     
//     LOGPRINTF( "@@@@@@@@@@ Exit: %f %f\n", vectornorm( M, ressigma ), vectornorm( pdim, resu  ) );
//     
//     
// }         



