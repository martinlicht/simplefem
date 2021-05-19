



#include "sparsesolver.hpp"


#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif






void ConjugateGradientSolverCSR( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float threshold,
    int print_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( residual );
    assert( threshold > 0 );
    assert( print_modulo >= -1 );
    
    Float* __restrict__ direction = new Float[N];
    Float* __restrict__ auxiliary = new Float[N];
    assert( direction );
    assert( auxiliary );
    
    Float r_r = notanumber;

    int K = 0;
    
    while( K < N ){
        
        bool restart_condition = ( K == 0 ); // or K % 1000 == 0;
        
        bool residual_seems_small = ( K != 0 ) and absolute(r_r) < threshold*threshold;

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
        
        /* printing information */

        if( print_modulo > 0 and K % print_modulo == 0 ) 
            printf("Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
                   K, N, (long double)(r_r), (long double) threshold*threshold );

        /* Check whether residual is small */
                
        bool residual_is_small = absolute(r_r) < threshold*threshold;
        
        if( residual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        // NOTE The calculation of d_r is reduced to r_r, which is already known.
        
        Float d_Ad = 0.;
        
        #pragma omp parallel for reduction(+:d_Ad) //d_r,
        for( int c = 0; c < N; c++ )
        {
            auxiliary[c] = 0.;
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                auxiliary[c] += csrvalues[ d ] * direction[ csrcolumns[d] ];
                    
            d_Ad += direction[c] * auxiliary[c];
        }
        
        
        bool denominator_is_unreasonable = not std::isfinite(d_Ad) or d_Ad < 0.;
        bool denominator_is_small    = sqrt(absolute(d_Ad)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) {
            if( print_modulo >= 0 ) printf( "Gradient energy is unreasonable with %.9Le\n", (long double)d_Ad );
            break;
        }
        
        if( denominator_is_small ) {
            if( print_modulo >= 0 ) printf( "Gradient energy is small with %.9Le while residual is %.9Le vs %.9Le\n", 
                                            (long double)d_Ad, (long double)r_r, (long double)threshold*threshold );
            break;
        }
        
        
        Float alpha = r_r / d_Ad;
    
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
        
        #pragma omp parallel 
        for( int c = 0; c < N; c++ )
            direction[c] = residual[c] + beta * direction[c];
        
        
        
        K++;
        
    }
    
    
    if( print_modulo >= 0 ) 
        printf( "Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
                K, N, (long double)(r_r), (long double) threshold*threshold );

    
    delete[] ( direction ); 
    delete[] ( auxiliary );

}



















void ConjugateGradientSolverCSR_DiagonalPreconditioner( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float threshold,
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
    assert( threshold > 0 );
    assert( print_modulo >= -1 );
    assert( precon );
    
    Float* __restrict__ direction = new Float[N];
    Float* __restrict__ zirconium = new Float[N];
    Float* __restrict__ auxiliary = new Float[N];
    assert( direction );
    assert( zirconium );
    assert( auxiliary );
    
    Float z_r = notanumber;

    int K = 0;
    
    while( K < N ){
        
        bool restart_condition = ( K == 0 ); // or K % 1000 == 0;
        
        bool preconresidual_seems_small = ( K != 0 ) and absolute(z_r) < threshold*threshold;

        if( restart_condition or preconresidual_seems_small ) {
            
            z_r = 0.;
        
            #pragma omp parallel for reduction(+:z_r)
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
            
        }
        
        /* printing information */

        if( print_modulo > 0 and K % print_modulo == 0 ) 
            printf("Residual after %d of max. %d iterations: %.9Le (%.9Le)\n",
                   K, N, (long double)(z_r), (long double) threshold*threshold );
        
        /* Check whether residual is small */
                
        bool preconresidual_is_small = absolute(z_r) < threshold*threshold;
        
        if( preconresidual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        // NOTE The calculation of d_r is reduced to r_r, which is already known.
        
        Float d_Ad = 0.;
        
        #pragma omp parallel for reduction(+:d_Ad) //d_r,
        for( int c = 0; c < N; c++ )
        {
            if( precon[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            auxiliary[c] = 0.;
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                auxiliary[c] += csrvalues[ d ] * direction[ csrcolumns[d] ];
                    
            d_Ad += direction[c] * auxiliary[c];
        }
        
        
        bool denominator_is_unreasonable = not std::isfinite(d_Ad) or d_Ad < 0.;
        bool denominator_is_small    = sqrt(absolute(d_Ad)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) {
            if( print_modulo >= 0 ) printf( "Gradient energy is unreasonable with %.9Le\n", (long double)d_Ad );
            break;
        }
        
        if( denominator_is_small ) {
            if( print_modulo >= 0 ) printf( "Gradient energy is small with %.9Le while precon-residual is %.9Le vs %.9Le\n", 
                                            (long double)d_Ad, (long double)z_r, (long double)threshold*threshold );
            break;
        }
        
        
        Float alpha = z_r / d_Ad;
    
        Float z_r_new = 0.;
        
        #pragma omp parallel for reduction(+:z_r_new) //r_r_old
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
        
        #pragma omp parallel 
        for( int c = 0; c < N; c++ ) {
            
            if( precon[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            direction[c] = zirconium[c] + beta * direction[c];
            
        }
        
        K++;
        
    }
    
    
    if( print_modulo >= 0 ) 
        printf( "Residual after %d of max. %d iterations: %.9Le (%.9Le)\n",
                K, N, (long double)(z_r), (long double) threshold*threshold );

    
    delete[] ( direction ); 
    delete[] ( zirconium ); 
    delete[] ( auxiliary );

}





























void ConjugateGradientSolverCSR_SSOR( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ residual,
    const Float threshold,
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
    assert( threshold > 0 );
    assert( print_modulo >= -1 );
    assert( diagonal );
    
    Float* __restrict__ direction = new Float[N];
    Float* __restrict__ zirconium = new Float[N];
    Float* __restrict__ auxiliary = new Float[N];
    Float* __restrict__ mittlerer = new Float[N];
    assert( direction );
    assert( zirconium );
    assert( mittlerer );
    assert( auxiliary );
    
    for( int i = 0; i < N; i++ )
        direction[i] = zirconium[i] = auxiliary[i] = mittlerer[i] = 0.;
    
    Float z_r = notanumber;

    int K = 0;
    
    while( K < N ){
        
        bool restart_condition = ( K == 0 ); // or K % 1000 == 0;
        
        bool preconresidual_seems_small = false and absolute(z_r) < threshold*threshold;

        if( restart_condition or preconresidual_seems_small ) {
            
            for( int c = 0; c < N; c++ ) {
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
                residual[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    residual[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
            }
            
            // inv( L^t + D/omega ) * D * inv( L + D/omega )
            
            for( int c = 0; c < N; c++ ) {
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                Float aux = residual[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    if( csrcolumns[d] < c )
                        aux -= csrvalues[ d ] * mittlerer[ csrcolumns[d] ];
                
                mittlerer[c] = aux * omega /  diagonal[c];
                
            }
        
            for( int c = 0; c < N; c++ ) {
            
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
                
                mittlerer[c] *= diagonal[c];
                
                mittlerer[c] *= ( 2. - omega ) / omega;
            }
            
            for( int c = N-1; c >= 0; c-- ) {
                
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
                Float aux = mittlerer[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    if( csrcolumns[d] > c )
                        aux -= csrvalues[ d ] * zirconium[ csrcolumns[d] ];
                
                zirconium[c] = aux * omega /  diagonal[c];
            
            }
            
            z_r = 0.;
        
            for( int c = 0; c < N; c++ ) {
                            
                if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
                direction[c] = zirconium[c];
                       
                z_r += zirconium[c] * residual[c];
            
            }
            
        }
        
        /* printing information */

        if( print_modulo > 0 and K % print_modulo == 0 ) 
            printf( "Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
                    K, N, (long double)(z_r), (long double) threshold*threshold );
        
        /* Check whether residual is small */
                
        bool preconresidual_is_small = absolute(z_r) < threshold*threshold;
        
        if( preconresidual_is_small )
            break;


        /* now the main work of the entire algorithm */
        
        // NOTE The calculation of d_r is reduced to r_r, which is already known.
        
        Float d_Ad = 0.;
        
        for( int c = 0; c < N; c++ )
        {
            auxiliary[c] = 0.;
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                auxiliary[c] += csrvalues[ d ] * direction[ csrcolumns[d] ];
                    
            d_Ad += direction[c] * auxiliary[c];
        }
        
        
        bool denominator_is_unreasonable = not std::isfinite(d_Ad) or d_Ad < 0.;
        bool denominator_is_small    = sqrt(absolute(d_Ad)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) {
            if( print_modulo >= 0 ) printf( "Gradient energy is unreasonable with %.9Le\n", (long double)d_Ad );
            break;
        }
        
        if( denominator_is_small ) {
            if( print_modulo >= 0 ) printf( "Gradient energy is small with %.9Le while precon-residual is %.9Le vs %.9Le\n", 
                                            (long double)d_Ad, (long double)z_r, (long double)threshold*threshold );
            break;
        }
        
        
        Float alpha = z_r / d_Ad;
    
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
        
        Float beta = z_r_new / z_r;
        
        z_r = z_r_new;
        
        for( int c = 0; c < N; c++ ) {
            
            if( diagonal[c] == 0. ) continue; // NOTE: guard against shadowed variables 
            
            direction[c] = zirconium[c] + beta * direction[c];
            
        }
        
        K++;
        
    }
    
    if( print_modulo >= 0 ) 
        printf( "Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
                K, N, (long double)(z_r), (long double) threshold*threshold );

    
    delete[] ( direction ); 
    delete[] ( zirconium ); 
    delete[] ( mittlerer ); 
    delete[] ( auxiliary );

}









































void ConjugateResidualSolverCSR( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ res,
    const Float threshold,
    int print_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( res );
    assert( threshold > 0 );
    assert( print_modulo >= -1 );
    
    Float* __restrict__  dir = new Float[N];
    Float* __restrict__ Adir = new Float[N];
    Float* __restrict__ Ares = new Float[N];
    assert(  dir );
    assert( Adir );
    assert( Ares );
    
    Float* __restrict__  vil = new Float[N];
    assert( vil );
    
    
    Float Ad_r  = notanumber;
    Float Ad_Ad = notanumber;
    
    int K = 0;
    
    while( K < N ){
        
        bool restart_condition = ( K == 0 ); // or K % 1000 == 0;
        
        bool residualenergy_seems_small = ( K != 0 ) and absolute(Ad_r) < threshold*threshold;

        if( restart_condition or residualenergy_seems_small ) {
            
            #pragma omp parallel for
            for( int c = 0; c < N; c++ ) {
                
                res[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    res[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                dir[c] = res[c]; // this line seems to slow down performance ....
                
            }
            
            Ad_r  = 0.;
            Ad_Ad = 0.;
            
            #pragma omp parallel for reduction(+:Ad_r,Ad_Ad)
            for( int c = 0; c < N; c++ ) {
                
                Adir[c] = 0.;
                
                for( int e = csrrows[c]; e < csrrows[c+1]; e++ )
                    Adir[c] += csrvalues[ e ] * dir[ csrcolumns[e] ];
                
                Ares[c] = Adir[c];
                
                Ad_r  += Adir[c] *  res[c];
                Ad_Ad += Adir[c] * Adir[c];
                
            }
            
            
        }
        
        /* printing information */

        if( print_modulo > 0 and K % print_modulo == 0 ) 
            printf( "Square Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
                    K, N, (long double)(Ad_r), (long double) threshold*threshold );

        
        /* Check whether res is small */
                
        bool residualenergy_is_small = absolute(Ad_r) < threshold*threshold;
        
        if( residualenergy_is_small )
            break;
        
        
        
        bool denominator_is_unreasonable = not std::isfinite(Ad_Ad) or Ad_Ad < 0.;
        bool denominator_is_small    = sqrt(absolute(Ad_Ad)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) {
            if( print_modulo >= 0 ) printf( "Gradient double energy is unreasonable with %.9Le\n", (long double)Ad_Ad );
            break;
        }
        
        if( denominator_is_small ) {
            if( print_modulo >= 0 ) printf( "Gradient double energy is small with %.9Le while precon-residual is %.9Le vs %.9Le\n", 
                                    (long double)Ad_Ad, (long double)Ad_r, (long double)threshold*threshold );
            break;
        }
        
        
        

        /* now the main work of the entire algorithm */
        
        Float alpha = Ad_r / Ad_Ad;
        
        Float new_Ar_r = 0.;
        
        #pragma omp parallel for reduction(+:new_Ar_r)
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
        #pragma omp parallel for reduction(+:Ad_Ad)
        for( int c = 0; c < N; c++ )
        {
            
            dir[c]  =  res[c] + beta *  dir[c];
            
            Adir[c] = Ares[c] + beta * Adir[c];
            
            Ad_Ad += Adir[c] * Adir[c];
        }
        
        Ad_r = new_Ar_r;
                
        K++;
        
    }
    
    if( print_modulo >= 0 ) 
        printf( "Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
                K, N, (long double)(Ad_r), (long double) threshold*threshold );

    
    delete[] (  dir );
    delete[] ( Adir );
    delete[] ( Ares );

    delete[] ( vil );

}




























void ConjugateResidualSolverCSR_textbook( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ res,
    const Float threshold,
    int print_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( res );
    assert( threshold > 0 );
    assert( print_modulo >= -1 );
    
    Float* __restrict__  dir = new Float[N];
    Float* __restrict__ Adir = new Float[N];
    Float* __restrict__ Ares = new Float[N];
    assert(  dir );
    assert( Adir );
    assert( Ares );
    
    Float* __restrict__  vil = new Float[N];
    assert( vil );
    
    
    Float Ar_r  = notanumber;
    Float Ad_Ad = notanumber;

    int K = 0;
    
    while( K < N ){
        
        bool restart_condition = ( K == 0 ); // or K % 1000 == 0;
        
        bool residualenergy_seems_small = ( K != 0 ) and absolute(Ar_r) < threshold*threshold;

        if( restart_condition or residualenergy_seems_small ) {
            
            #pragma omp parallel for
            for( int c = 0; c < N; c++ ) {
                
                res[c] = b[c];
                
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    res[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
                
                dir[c] = res[c]; // this line seems to slow down performance ....
                
            }
            
            Ar_r  = 0.;
            Ad_Ad = 0.;
            
            #pragma omp parallel for reduction(+:Ar_r,Ad_Ad)
            for( int c = 0; c < N; c++ ) {
                
                Adir[c] = 0.;
                
                for( int e = csrrows[c]; e < csrrows[c+1]; e++ )
                    Adir[c] += csrvalues[ e ] * dir[ csrcolumns[e] ];
                
                Ares[c] = Adir[c];
                
                Ar_r  += Ares[c] *  res[c];
                Ad_Ad += Adir[c] * Adir[c];
                
            }
            
            
        }
        
        /* printing information */

        if( print_modulo > 0 and K % print_modulo == 0 )
            printf( "Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
                    K, N, (long double)(Ar_r), (long double) threshold*threshold );
        
        /* Check whether res is small */
                
        bool residualenergy_is_small = absolute(Ar_r) < threshold*threshold;
        
        if( residualenergy_is_small )
            break;


        bool denominator_is_unreasonable = not std::isfinite(Ad_Ad) or Ad_Ad < 0.;
        bool denominator_is_small    = sqrt(absolute(Ad_Ad)) < machine_epsilon;
        
        if( denominator_is_unreasonable ) {
            if( print_modulo >= 0 ) printf( "Gradient energy is unreasonable with %.9Le\n", (long double)Ar_r );
            break;
        }
        
        if( denominator_is_small ) {
            if( print_modulo >= 0 ) printf( "Gradient energy is small with %.9Le while precon-residual is %.9Le vs %.9Le\n", 
                                            (long double)Ad_Ad, (long double)Ar_r, (long double)threshold*threshold );
            break;
        }
        
        
        /* now the main work of the entire algorithm */
        
        Float alpha = Ar_r / Ad_Ad;
        
        Float new_Ar_r = 0.;
        
        #pragma omp parallel for reduction(+:new_Ar_r)
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
        #pragma omp parallel for reduction(+:Ad_Ad)
        for( int c = 0; c < N; c++ )
        {
            
            dir[c]  =  res[c] + beta *  dir[c];
            
            Adir[c] = Ares[c] + beta * Adir[c];
            
            Ad_Ad += Adir[c] * Adir[c];
        }
        
        Ar_r = new_Ar_r;
                
        K++;
        
    }
    
    if( print_modulo >= 0 ) 
        printf( "Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", 
                K, N, (long double)(Ar_r), (long double) threshold*threshold );

    
    delete[] (  dir );
    delete[] ( Adir );
    delete[] ( Ares );

    delete[] ( vil );

}




































                      




void MINRESCSR( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ res,
    const Float threshold,
    int print_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( res );
    assert( threshold > 0 );
    assert( print_modulo >= -1 );
    
    Float* __restrict__ v0 = new Float[N];
    Float* __restrict__ v1 = new Float[N];
    Float* __restrict__ w0 = new Float[N];
    Float* __restrict__ w1 = new Float[N];
    
    assert( v0 );
    assert( v1 );
    assert( w0 );
    assert( w1 );
    
    Float* __restrict__ vn = new Float[N];
    Float* __restrict__ wn = new Float[N];
    Float* __restrict__  p = new Float[N];
    
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

    while( K < N ){
        
        
        bool restart_condition = (K == 0);
        
        bool residual_seems_small = ( K != 0 ) and (eta < threshold);
        
        if( restart_condition or residual_seems_small ) {
            
            Float gamma_sq = 0.;
            
            #pragma omp parallel for reduction(+:gamma_sq)
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
            
            #pragma omp parallel for 
            for( int c = 0; c < N; c++ ) 
                v1[c] /= gamma;
            
            s0 = s1 = 0.;
            c0 = c1 = 1.;
            
            eta = gamma;
            
        }
        
        bool residual_is_small = (eta < threshold);
        
        if( residual_is_small )
            break;

        Float temp = gamma;
//         if( K % 100 == 0 )
//         LOG << K << space << temp << space << eta << space << temp/eta;
            
        {
            
            
            Float delta = 0.;
            
            #pragma omp parallel for reduction(+:delta)
            for( int c = 0; c < N; c++ ) {
                
                p[c] = 0.;
                for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
                    p[c] += csrvalues[ d ] * v1[ csrcolumns[d] ];
                
                delta += p[c] * v1[c];
                
            }
            
            assert( delta > 0. );
            
            Float gamma_n_sq = 0.;
            #pragma omp parallel for reduction(+:gamma_n_sq)
            for( int c = 0; c < N; c++ ) {
                
                vn[c] = p[c] - delta * v1[c] - gamma * v0[c];
                
                gamma_n_sq += vn[c] * vn[c];
                
            }
            
            assert( gamma_n_sq > 0. );
 
            Float gamma_n = std::sqrt(gamma_n_sq);
            
            #pragma omp parallel for 
            for( int c = 0; c < N; c++ )
                vn[c] /= gamma_n;
            
            Float alpha_0 = c1 * delta - c0 * s1 * gamma;
            Float alpha_1 = std::sqrt( alpha_0 * alpha_0 + gamma_n * gamma_n );
            Float alpha_2 = s1 * delta + c0 * c1 * gamma;
            Float alpha_3 = s0 * gamma;
 
            assert( alpha_1 > 0. );

            Float cn = alpha_0 / alpha_1;
            Float sn = gamma_n / alpha_1;
            
            #pragma omp parallel for
            for( int c = 0; c < N; c++ ) {
                wn[c] = ( v1[c] - alpha_2 * w1[c] - alpha_3 * w0[c] ) / alpha_1;
                x[c] = x[c] + cn * eta * wn[c];
            }
            
            eta = - sn * eta;
            
//             LOG << "\t" << alpha_0 << space << alpha_1 << space << alpha_2 << space << alpha_3;
//             LOG << "\t" << cn << space << sn << space << eta;
            
            
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

        
        if( print_modulo > 0 and K % print_modulo == 0 )
            printf("Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", K, N, (long double)eta, (long double) threshold );
        
        K++;
        
    }
    
    if( print_modulo >= 0 ) 
        printf("Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", K, N, (long double)eta, (long double) threshold );

    
    delete[] ( v0 );
    delete[] ( v1 );
    delete[] ( w0 );
    delete[] ( w1 );
    
    delete[] ( vn );
    delete[] ( wn );
    delete[] (  p );
    
}








// This one is taken from Wikipedia, I don't know what it is.

void WHATEVER( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ res,
    const Float threshold,
    int print_modulo
) {
    
    assert( N > 0 );
    assert( x );
    assert( b );
    assert( csrrows );
    assert( csrcolumns );
    assert( csrvalues );
    assert( res );
    assert( threshold > 0 );
    assert( print_modulo >= -1 );
    
    Float* __restrict__  r = new Float[N];
    Float* __restrict__ p0 = new Float[N];
    Float* __restrict__ p1 = new Float[N];
    Float* __restrict__ p2 = new Float[N];
    Float* __restrict__ s0 = new Float[N];
    Float* __restrict__ s1 = new Float[N];
    Float* __restrict__ s2 = new Float[N];
    
    assert(  r );
    assert( p0 );
    assert( p1 );
    assert( p2 );
    assert( s0 );
    assert( s1 );
    assert( s2 );

    Float r_r = 0.;
    
    #pragma omp parallel for reduction(+:r_r)
    for( int c = 0; c < N; c++ ) {
        
        r[c] = b[c];
        
        for( int d = csrrows[c]; d < csrrows[c+1]; d++ )
            r[c] -= csrvalues[ d ] * x[ csrcolumns[d] ];
        
        p0[c] =  r[c];
        p1[c] = p0[c];
        
        r_r += r[c] * r[c];
                
    }
    
    #pragma omp parallel for
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
        
        #pragma omp parallel for reduction(+:r_s1,s0_s0)
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
        
        #pragma omp parallel for reduction(+:r_r)
        for( int c = 0; c < N; c++ ) {
            
            x[c] += alpha * p1[c];
            r[c] -= alpha * s1[c];
            
            r_r += r[c] * r[c];
        }
        
        if( std::sqrt(r_r) < threshold )
            break;
        
        Float s0_s1 = 0.;
        Float s1_s1 = 0.;
        
        #pragma omp parallel for reduction(+:s0_s1,s1_s1)
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
        
        #pragma omp parallel for reduction(+:s0_s2,s2_s2)
        for( int c = 0; c < N; c++ ) {
            
            p0[c] -= beta * p1[c];
            s0[c] -= beta * s1[c];
            
            s0_s2 += s0[c] * s2[c];
            s2_s2 += s2[c] * s2[c];
        }
        
        if( K > 0 )
        {
            
            Float gamma = s0_s2 / s2_s2;
            
            #pragma omp parallel for 
            for( int c = 0; c < N; c++ ) {
                
                p0[c] -= gamma * p2[c];
                s0[c] -= gamma * s2[c];
                
            }
            
        }
            
            
        
        if( print_modulo > 0 and K % print_modulo == 0 ) 
            printf("Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", K, N, (long double)std::sqrt(r_r), (long double) threshold );
        
        K++;
        
    }
    
    if( print_modulo >= 0 ) 
        printf("Residual after %d of max. %d iterations: %.9Le (%.9Le)\n", K, N, (long double)std::sqrt(r_r), (long double) threshold );

    
    delete[] (  r );
    delete[] ( p0 );
    delete[] ( p1 );
    delete[] ( p2 );
    delete[] ( s0 );
    delete[] ( s1 );
    delete[] ( s2 );
    
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
//     const Float threshold = threshold;
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
//         if( residualnorm < threshold ) {
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
//     const Float threshold = threshold;
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
//         if( residualnorm < threshold ) {
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
//     printf( "@@@@@@@@@@ Exit: %f %f\n", vectornorm( M, ressigma ), vectornorm( pdim, resu  ) );
//     
//     
// }         



