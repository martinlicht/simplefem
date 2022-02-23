
#include <utility>

#include "iterativesolver.hpp"

#include "../operators/floatvector.hpp"


inline const bool cpp_restart_on_full_dimension = false;
inline const bool cpp_restart_before_finish     = false;










ConjugateGradientMethod::ConjugateGradientMethod( const LinearOperator& op )
: IterativeSolver(), A( op )
{
    this->max_iteration_count = op.getdimin();
    ConjugateGradientMethod::check();
}

ConjugateGradientMethod::~ConjugateGradientMethod()
{
    ConjugateGradientMethod::check();
}

void ConjugateGradientMethod::check() const
{
    IterativeSolver::check();
    assert( A.getdimin() == A.getdimout() );
}

std::string ConjugateGradientMethod::text() const
{
    return "Solver: Conjugate Gradient Method";
}




  

  
void ConjugateGradientMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
    const int dimension = A.getdimin();

    /* Build up data */
    
    FloatVector  r( dimension, 0. );
    FloatVector  d( dimension, 0. );
    
    // avoid repeated allocation of these temporary vectors 
    FloatVector Ad( dimension, 0. );
    
    recent_iteration_count = 0;

    Float sigma_min_sq = b * ( A * b ) / (b*b);
    
    if( verbosity >= VerbosityLevel::verbose ) LOG << "START Conjugate Gradient iteration" << nl;
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and absolute( r * r ) < threshold*threshold;
        
        /* Start / Restart CRM process */
        if( restart_condition or ( cpp_restart_before_finish and residual_seems_small ) ) {
        
            r = b - A * x;
            d = r;

            if( verbosity >= VerbosityLevel::verbose )
                LOGPRINTF( "RESTARTED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(r*r), (long double)threshold );
            
        }

        bool residual_is_small = absolute( r * r ) < threshold*threshold;
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( verbosity >= VerbosityLevel::verbose and print_condition )
            LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(r*r), (long double)threshold );
        
        /* If exit condition met, exit */
        
        if( residual_is_small ) 
            break;
            
        /* Perform iteration step */
        {
            
            Ad = A * d;
            
            Float rr_old = r * r;
            Float Ad_d  = Ad * d;
            Float alpha = rr_old / Ad_d;

            assert( Ad_d >= 0 );

            // sigma_min_sq = std::min( sigma_min_sq, Ad_d / (d*d) );
            // LOG << "@" << recent_iteration_count << " : " << std::sqrt( rr_old / ( sigma_min_sq * (x*x) ) ) << " with eigenvalue bound " << std::sqrt(sigma_min_sq) << nl;

            bool denominator_is_unreasonable = not std::isfinite(Ad_d) or Ad_d < 0.;
            bool denominator_is_small    = sqrt(absolute(Ad_d)) < machine_epsilon;
            
            if( denominator_is_unreasonable ) {
                LOGPRINTF( "BREAKDOWN: Gradient energy is unreasonable with %.9Le\n", (long double)Ad_d );
                break;
            }
            
            if( denominator_is_small ) {
                LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(rr_old), (long double)threshold );
                LOGPRINTF( "WARNING: Gradient energy is small with %.9Le\n", (long double)Ad_d );
                break;
            }
            
            x = x + alpha * d;
            r = r - alpha * Ad;
            
            Float rr_new = r * r;
            Float beta = rr_new / rr_old;
            d = r + beta * d;
        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
    
    /* HOW DID WE FINISH ? */
    recent_deviation = absolute(r * r);
    
    if( verbosity >= VerbosityLevel::resultonly and print_modulo >= 0 ) 
        LOGPRINTF( "FINISHED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(recent_deviation), (long double)threshold ); 
        
    
}

  
  
/* 
 * 
 * Input: A, b, x
 *****************
 * r := b - A x
 * d := r
 * Ar = A r
 * Ad = A d
 * rho = r . Ar
 * loop until exit condition
 *   p = A * Ad
 *   alpha = rho / d . p
 *
 *   x += alpha d
 *   r -= alpha Ad
 *   Ar -= alpha p
 *   rho' = r. Ar
 *   beta = rho' / rho
 *   
 *   rho = rho'
 *   d = r + beta d
 *   Ad = Ar + beta Ad
 * 
 */  
  
  
  



/* 
 * 
 * Input: A, b, x
 *****************
 * r := b - A x
 * d := r
 * Ar = A r
 * Ad = A d
 * rho = r . Ar
 * loop until exit condition
 *   p = A * Ad
 *   alpha = rho / Ad . Ad
 *
 *   x += alpha d
 *   r -= alpha Ad
 *   Ar -= alpha p
 *   rho' = r. Ar
 *   beta = rho' / rho
 *   
 *   rho = rho'
 *   d = r + beta d
 *   Ad = Ar + beta Ad
 * 
 */  

 /* 
 * 
 * Input: A, b, x
 *****************
 * r := b - A x
 * d := r
 * Ar = A r
 * Ad = A d
 * rho = r . Ar
 * p = A * Ad 
 * alpha = rho / Ad . Ad
 * loop until exit condition
 *
 *   x += alpha d
 *   r -= alpha Ad
 *   Ar -= alpha p
 *   rho' = r. Ar
 *   beta = rho' / rho
 *   
 *   rho = rho'
 *   d = r + beta d
 *   Ad = Ar + beta Ad
 *   p = A * Ad
 *   alpha = rho / Ad . Ad
 * 
 */  
 


















ConjugateResidualMethod::ConjugateResidualMethod( const LinearOperator& op )
: IterativeSolver(), A( op )
{
    this->max_iteration_count = op.getdimin();
    ConjugateResidualMethod::check();
}

ConjugateResidualMethod::~ConjugateResidualMethod()
{
    ConjugateResidualMethod::check();
}

void ConjugateResidualMethod::check() const
{
    IterativeSolver::check();
    assert( A.getdimin() == A.getdimout() );
}

std::string ConjugateResidualMethod::text() const
{
    return "Solver: Conjugate Residual Method";
}



















void ConjugateResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    solve_robust( x, b );
}




  
void ConjugateResidualMethod::solve_explicit( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
    const int dimension = A.getdimin();

    /* Build up data */
    
    Float rAr = notanumber;

    FloatVector  r( dimension, 0. );
    FloatVector  d( dimension, 0. );
    
    FloatVector Ad( dimension, 0. );
    FloatVector Ar( dimension, 0. );
    
    FloatVector AAd( dimension, 0. );
    
    recent_iteration_count = 0;
    
    if( verbosity >= VerbosityLevel::verbose ) LOG << "START Conjugate Residual iteration" << nl;
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        /* Start / Restart CRM process */
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and ( absolute( r * r ) < threshold*threshold or absolute( rAr ) < threshold*threshold );
        
        if( restart_condition or ( cpp_restart_before_finish and residual_seems_small ) ) {
        
            /* r = b - A x */
            r = b - A * x;

            /* d = r */
            d = r; //d.copydatafrom( r );

            /* Ar = A r */
            Ar = A * r; // A.apply( Ar, (const FloatVector&) r );

            /* Ad = A d */
            Ad = A * d; // A.apply( Ad, (const FloatVector&) d );

            /* rho is r.A.r */
            rAr = Ar * r;
            if( rAr < 0. ) LOG << rAr << nl;
            assert( rAr >= 0. );

            if( verbosity >= VerbosityLevel::verbose )
                LOGPRINTF( "RESTARTED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(rAr), (long double)threshold );
            

            
        }

        bool residual_is_small = absolute( r * r ) < threshold*threshold or absolute( rAr ) < threshold*threshold; 
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( verbosity >= VerbosityLevel::verbose and print_condition )
            LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(rAr), (long double)threshold ); 
        
        /* If exit condition met, exit */
        if( residual_is_small ) 
            break;
            
        
        /* Perform iteration step */
        {

            if( rAr < 0. ) LOG << rAr << nl;
            assert( rAr >= 0. );
            
            AAd = A * Ad;

            /*  alpha = r.A.r / d.AAd */
            Float Ad_Ad = Ad * Ad;
            Float alpha = rAr / Ad_Ad;

            bool denominator_is_unreasonable = not std::isfinite(Ad_Ad) or Ad_Ad < 0.;
            bool denominator_is_small    = sqrt(absolute(Ad_Ad)) < machine_epsilon;
            
            if( denominator_is_unreasonable ) {
                LOGPRINTF( "BREAKDOWN: Gradient double energy is unreasonable with %.9Le\n", (long double)Ad_Ad );
                break;
            }
            
            if( denominator_is_small ) {
                LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(rAr), (long double)threshold ); 
                LOGPRINTF( "WARNING: Gradient double energy is small with %.9Le\n", (long double)Ad_Ad );
                break;
            }
            
                        
            x += alpha * d;

            r -= alpha * Ad; //= b - A * x;

            Ar -= alpha * AAd;

            Float tau = rAr;
            rAr = r * Ar; //r * ( A * r ); //r * Ar;
            Float beta = rAr / tau;
            
            if( rAr < 0. ) {
                LOGPRINTF( "BREAKDOWN: Residual energy is unreasonable with %.9Le\n", (long double)rAr );
                rAr = 0.; // TODO to avoid useless bug 
                break;
            }

            assert( rAr >= 0. );
            
            d = r + beta * d;

            Ad = Ar + beta * Ad;
            
        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
    
    /* HOW DID WE FINISH ? */
    recent_deviation = rAr;
    
    if( verbosity >= VerbosityLevel::resultonly and print_modulo >= 0 ) 
        LOGPRINTF( "FINISHED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(recent_deviation), (long double)threshold ); 

    
}
  

  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
void ConjugateResidualMethod::solve_robust( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
    const int dimension = A.getdimin();

    FloatVector  r( dimension, 0. );
    FloatVector  d( dimension, 0. );
    FloatVector Ar( dimension, 0. );
    FloatVector Ad( dimension, 0. );
    FloatVector  p( dimension, 0. );
    
    recent_iteration_count = 0;
    
    if( verbosity >= VerbosityLevel::verbose ) LOG << "START Conjugate Residual iteration" << nl;
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and ( absolute( r * r ) < threshold*threshold or absolute( r * Ar ) < threshold*threshold );
        // first criterion is not in fast 

        if( restart_condition or ( cpp_restart_before_finish and residual_seems_small ) ) {
        
            r  = b - A * x;
            d  = r;
            
            Ar = A * r;
            Ad = Ar;

            // fast and explicit: Ar_r = Ar * r;

            if( verbosity >= VerbosityLevel::verbose )
                LOGPRINTF( "RESTARTED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt( r * Ar ), (long double)threshold );

        }

        bool residual_is_small = absolute( r * r ) < threshold*threshold or absolute( r * Ar ) < threshold*threshold; 
        // first criterion is in explicit but not in fast (for speed reasons)
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( verbosity >= VerbosityLevel::verbose and print_condition )
            LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt( Ar * r ), (long double)threshold ); 
        
        /* If exit condition met, exit */
        
        if( residual_is_small ) 
            break;
            
        {
            p = A * Ad;
            Float Ad_Ad = Ad * Ad; assert( Ad_Ad >= 0 ); 
            
            Float Ad_r  = Ad * r; // fast/explicit uses and maintains Ar_r
            Float alpha = Ad_r / Ad_Ad;
            
            bool denominator_is_unreasonable = not std::isfinite(Ad_Ad) or Ad_Ad < 0.;
            bool denominator_is_small        = sqrt(absolute(Ad_Ad)) < machine_epsilon;
            
            if( denominator_is_unreasonable ) {
                LOGPRINTF( "BREAKDOWN: Gradient double energy is unreasonable with %.9Le\n", (long double)Ad_Ad );
                break;
            }
            
            if( denominator_is_small ) {
                LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt( Ar * r ), (long double)threshold ); 
                LOGPRINTF( "WARNING: Gradient double energy is small with %.9Le\n", (long double)Ad_Ad );
                break;
            }
            
            x  =  x + alpha * d;
            r  =  r - alpha * Ad;
            Ar = Ar - alpha * p;
            
            Float Ar_r_new = Ar * r;
            Float beta = Ar_r_new / Ad_r; // fast uses and maintains Ar_r // TODO: abort if Ad_r is negative?
            
            // assert( Ar_r_new >= 0. );

            d  =  r + beta *  d;
            Ad = Ar + beta * Ad;
            
        }
        
        recent_iteration_count++;
    }
    
    
    recent_deviation = absolute( r * Ar );
    
    if( verbosity >= VerbosityLevel::resultonly and print_modulo >= 0 ) 
        LOGPRINTF( "FINISHED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(recent_deviation), (long double)threshold ); 
    
}






















void ConjugateResidualMethod::solve_fast( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
    const int dimension = A.getdimin();

    FloatVector  r( dimension, 0. );
    FloatVector  d( dimension, 0. );
    FloatVector Ar( dimension, 0. );
    FloatVector Ad( dimension, 0. );
    FloatVector  p( dimension, 0. );
    
    recent_iteration_count = 0;
    
    Float Ar_r = notanumber;
    
    if( verbosity >= VerbosityLevel::verbose ) LOG << "START Conjugate Residual iteration" << nl;
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and absolute( Ar * r ) < threshold*threshold;
        
        if( restart_condition or ( cpp_restart_before_finish and residual_seems_small ) ) {
        
            r  = b - A * x;
            d  = r;
            
            Ar = A * r;
            Ad = Ar;
            
            Ar_r = Ar * r;

            if( verbosity >= VerbosityLevel::verbose )
                LOGPRINTF( "RESTARTED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(Ar_r), (long double)threshold );

        }

        bool residual_is_small = absolute( Ar * r ) < threshold*threshold; 
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( verbosity >= VerbosityLevel::verbose and print_condition )
            LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt( Ar * r ), (long double)threshold );
        
        /* If exit condition met, exit */
        
        if( residual_is_small ) 
            break;
            
        {
            
            p = A * Ad;
            Float Ad_Ad = Ad * Ad; assert( Ad_Ad >= 0 ); 
            
            Float alpha = Ar_r / Ad_Ad;
            
            bool denominator_is_unreasonable = not std::isfinite(Ad_Ad) or Ad_Ad < 0.;
            bool denominator_is_small    = sqrt(absolute(Ad_Ad)) < machine_epsilon;
            
            if( denominator_is_unreasonable ) {
                LOGPRINTF( "BREAKDOWN: Gradient double energy is unreasonable with %.9Le\n", (long double)Ad_Ad );
                break;
            }
            
            if( denominator_is_small ) {
                LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt( Ar * r ), (long double)threshold );
                LOGPRINTF( "WARNING: Gradient double energy is small with %.9Le\n", (long double)Ad_Ad );
                break;
            }
            
            x  =  x + alpha * d;
            r  =  r - alpha * Ad;
            Ar = Ar - alpha * p;
            
            Float Ar_r_new = Ar * r;
            Float beta = Ar_r_new / Ar_r;
            
            // assert( Ar_r_new >= 0. );

            d  =  r + beta *  d;
            Ad = Ar + beta * Ad;
            
            Ar_r = Ar_r_new;
            
        }
        
        recent_iteration_count++;
    }
    
    
    recent_deviation = absolute( Ar * r );
    
    if( verbosity >= VerbosityLevel::resultonly and print_modulo >= 0 ) 
        LOGPRINTF( "FINISHED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(recent_deviation), (long double)threshold );
    
}


  
/* 
 * 
 * Input: A, b, x
 *****************
 * r := b - A x
 * d := r
 * Ar = A r
 * Ad = A d
 * rho = r . Ar
 * loop until exit condition
 *   p = A * Ad
 *   alpha = rho / d . p
 *
 *   x += alpha d
 *   r -= alpha Ad
 *   Ar -= alpha p
 *   rho' = r. Ar
 *   beta = rho' / rho
 *   
 *   rho = rho'
 *   d = r + beta d
 *   Ad = Ar + beta Ad
 * 
 */  
  
  
  



/* 
 * 
 * Input: A, b, x
 *****************
 * r := b - A x
 * d := r
 * Ar = A r
 * Ad = A d
 * rho = r . Ar
 * loop until exit condition
 *   p = A * Ad
 *   alpha = rho / Ad . Ad
 *
 *   x += alpha d
 *   r -= alpha Ad
 *   Ar -= alpha p
 *   rho' = r. Ar
 *   beta = rho' / rho
 *   
 *   rho = rho'
 *   d = r + beta d
 *   Ad = Ar + beta Ad
 * 
 */  

 /* 
 * 
 * Input: A, b, x
 *****************
 * r := b - A x
 * d := r
 * Ar = A r
 * Ad = A d
 * rho = r . Ar
 * p = A * Ad 
 * alpha = rho / Ad . Ad
 * loop until exit condition
 *
 *   x += alpha d
 *   r -= alpha Ad
 *   Ar -= alpha p
 *   rho' = r. Ar
 *   beta = rho' / rho
 *   
 *   rho = rho'
 *   d = r + beta d
 *   Ad = Ar + beta Ad
 *   p = A * Ad
 *   alpha = rho / Ad . Ad
 * 
 */  

 
 

























  
PreconditionedConjugateResidualMethod::PreconditionedConjugateResidualMethod( const LinearOperator& op, const LinearOperator& M )
: IterativeSolver(), A( op ), M( M )
{
    this->max_iteration_count = op.getdimin();
    PreconditionedConjugateResidualMethod::check();
}

PreconditionedConjugateResidualMethod::~PreconditionedConjugateResidualMethod()
{
    PreconditionedConjugateResidualMethod::check();
}
  
void PreconditionedConjugateResidualMethod::check() const
{
    IterativeSolver::check();
    assert( A.getdimin() == A.getdimout() );
    assert( M.getdimin() == M.getdimout() );
    assert( A.getdimin() == M.getdimout() );
}

std::string PreconditionedConjugateResidualMethod::text() const
{
    return "Solver: Preconditioned Conjugate Residual Method.";
}





void PreconditionedConjugateResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );

    const int dimension = A.getdimin();

    /* Build up data */
    
    Float rMAMr = notanumber;
    
    FloatVector   r( dimension, 0. );
    FloatVector   p( dimension, 0. ); 
    
    FloatVector  Mr( dimension, 0. );
    FloatVector  Mp( dimension, 0. );
    FloatVector AMr( dimension, 0. );
    FloatVector AMp( dimension, 0. );
        
    // avoid repeated allocation of these temporary vectors 
    FloatVector  MAMp( dimension, 0. );    
    FloatVector AMAMp( dimension, 0. );
        
    recent_iteration_count = 0;
    
    if( verbosity >= VerbosityLevel::verbose ) LOG << "START Preconditioned Conjugate Residual iteration" << nl;
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and absolute( rMAMr ) < threshold*threshold;
        
        /* Start / Restart PCRM process */
        if( restart_condition or ( cpp_restart_before_finish and residual_seems_small ) ) {
        
            {
                
                /* x = initial guess */
                
                /* r = b - A x */
                /* p = r */
                
                r = b - A * x;
                p = M * ( A * r );
                
                /* Mr = M r */
                /* Mp = M p */
                Mr = M * r;
                Mp = M * p;

                /* Ar = A r */
                /* Ap = A p */
                AMr = A * Mr;
                AMp = A * Mp;

                /* rho is Mr.A.Mr */
                rMAMr = Mr * AMr;
                
                if( verbosity >= VerbosityLevel::verbose )
                    LOGPRINTF( "RESTARTED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(rMAMr), (long double)threshold );

            }

        }

        /* If exit condition met, exit */

        bool residual_is_small = absolute( rMAMr ) < threshold*threshold;
        
        if( residual_is_small ) 
            break;
            
        /* Print information */
        
        if( print_modulo > 0 and recent_iteration_count % print_modulo == 0 ) 
            LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(rMAMr), (long double)threshold );

        /* Perform iteration step */
        {
            
            MAMp  = M *  AMp;
            AMAMp = A * MAMp;

            Float AMp_MAMp = AMp * MAMp;
            Float alpha = rMAMr / ( AMp_MAMp );

            bool denominator_is_unreasonable = not std::isfinite(AMp_MAMp) or AMp_MAMp < 0.;
            bool denominator_is_small    = sqrt(absolute(AMp_MAMp)) < machine_epsilon;
            
            if( denominator_is_unreasonable ) {
                LOGPRINTF( "BREAKDOWN: Gradient double energy is unreasonable with %.9Le\n", (long double)AMp_MAMp );
                break;
            }
            
            if( denominator_is_small ) {
                LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(Mr * AMr), (long double)threshold );
                LOGPRINTF( "WARNING: Gradient double energy is small with %.9Le\n", (long double)AMp_MAMp );
                break;
            }
            
            x += alpha * Mp;

            r  -= alpha *   AMp; // not necessary, see also below
            Mr -= alpha *  MAMp;
            AMr -= alpha * AMAMp;

            Float new_rMAMr = Mr * AMr;

            Float beta = new_rMAMr / rMAMr;

            p =   r + beta *   p; // as such, completely useless
            Mp =  Mr + beta *  Mp;
            AMp = AMr + beta * AMp;

            rMAMr = new_rMAMr;

        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
        

    /* HOW DID WE FINISH ? */

    recent_deviation = rMAMr;
    
    if( verbosity >= VerbosityLevel::resultonly and print_modulo >= 0 ) 
        LOGPRINTF( "FINISHED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(recent_deviation), (long double)threshold );

    
}
  








































MinimumResidualMethod::MinimumResidualMethod( const LinearOperator& op )
: IterativeSolver(), A( op )
{
    this->max_iteration_count = op.getdimin();
    MinimumResidualMethod::check();
}

MinimumResidualMethod::~MinimumResidualMethod()
{
    MinimumResidualMethod::check();
}

void MinimumResidualMethod::check() const
{
    IterativeSolver::check();
    assert( A.getdimin() == A.getdimout() );
}

std::string MinimumResidualMethod::text() const
{
    return "Solver: Minimum Residual Method";
}



  
void MinimumResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    assert( x.isfinite() and b.isfinite() );
    
    const int dimension = A.getdimin();

    /* Build up data */
    
    Float rr = notanumber;
    
    FloatVector  r( dimension, 0. );
    FloatVector p0( dimension, 0. );
    FloatVector s0( dimension, 0. );
    FloatVector p1( dimension, 0. );
    FloatVector s1( dimension, 0. );
    
    // auxiliary data to avoid reallocation
    FloatVector p2( dimension, 0. );
    FloatVector s2( dimension, 0. );
    
    /* Begin iteration */
    
    recent_iteration_count = 0;
    
//     Float recent_alpha = notanumber;
    
    if( verbosity >= VerbosityLevel::verbose ) LOG << "START Minimal Residual iteration" << nl;
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and absolute( rr ) < threshold*threshold;
        
        /* Start / Restart MinimumResidualMethod process */
        if( restart_condition or ( cpp_restart_before_finish and residual_seems_small ) ) {
        
            r = b - A * x;
            
            rr = r * r;
            
            assert( std::isfinite(rr) );
            if( rr < threshold*threshold )
                break;
            
            p0 = r;
            s0 = A * r;
            
            
            Float s0_s0 = s0 * s0;
            
            assert( std::isfinite(s0_s0) );
            if( s0_s0 < threshold*threshold )
                break;
            
            p0 /= std::sqrt(s0_s0);
            s0 /= std::sqrt(s0_s0);
            
            Float alpha0 = r * s0;
            
            x = x + alpha0 * p0;
            r = r - alpha0 * s0;
            
            rr = r * r;
            
            assert( std::isfinite(rr) );
            if( rr < threshold*threshold )
                break;
            
            
            p1 = r; 
            s1 = A * r;
            
            Float beta0 = s0 * s1;
            
            p1 = p1 - beta0 * p0;
            s1 = s1 - beta0 * s0;
            
            Float s1_s1 = s1 * s1;
            
            assert( std::isfinite(s1_s1) );
            if( s1_s1 < threshold*threshold )
                break;
            
            p1 /= std::sqrt(s1_s1);
            s1 /= std::sqrt(s1_s1);
            
            Float alpha1 = r * s1;
            
            x = x + alpha1 * p1;
            r = r - alpha1 * s1;
            
            rr = r * r;
            
//             recent_alpha = alpha1;

            if( verbosity >= VerbosityLevel::verbose )
                LOGPRINTF( "RESTARTED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(rr), (long double)threshold );

            
        }

        bool residual_is_small = absolute( rr ) < threshold*threshold; 
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( verbosity >= VerbosityLevel::verbose and print_condition )
            LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(rr), (long double)threshold );
        
        /* If exit condition met, exit */
        
        if( residual_is_small ) 
            break;
            
            
        /* Perform iteration step */
        {

            p2 = r;
            s2 = A * r;
            
            for( int i = 0; i < 2; i++ ){
                
                Float beta0 = s0 * s2;
                
                p2 = p2 - beta0 * p0;
                s2 = s2 - beta0 * s0;
                
//                 LOG << beta0 / ( A * r ).norm_sq() << nl;
                
                
                Float beta1 = s1 * s2;
                
                p2 = p2 - beta1 * p1;
                s2 = s2 - beta1 * s1;
                
            }
            
            
            Float s2_s2 = s2 * s2;
            
            assert( std::isfinite(s2_s2) ); 
            if( s2_s2 < machine_epsilon ) {
                LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(rr), (long double) threshold );
                LOGPRINTF( "WARNING: Norm of search direction is small with %.9Le\n", (long double) sqrt(s2_s2) );
                LOGPRINTF( "WARNING: Machine Epsilon: %.9Le\n", (long double) machine_epsilon );
                // break;
            }
            
            Float mysqrt = std::sqrt(s2_s2);
            p2 /= mysqrt;
            s2 /= mysqrt;
            
            Float alpha2 = r * s2;
            
            x = x + alpha2 * p2;
            r = r - alpha2 * s2;
            
//             r = b - A * x;
            
            rr = r * r;
            
//             recent_alpha = alpha2;
            
            assert( std::isfinite(rr) );
            if( sqrt(rr) < threshold )
                break;
            
            
            std::swap( p0, p1 ); std::swap( p1, p2 );
            std::swap( s0, s1 ); std::swap( s1, s2 );
//             p0 = p1; p1 = p2;
//             s0 = s1; s1 = s2;
            
        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
    
    /* HOW DID WE FINISH ? */
    recent_deviation = absolute(rr);
    
    if( verbosity >= VerbosityLevel::resultonly and print_modulo >= 0 ) 
        LOGPRINTF( "FINISHED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(recent_deviation), (long double)threshold );
    
}











/*  
function [x,r] = minres(A,b,x0,maxit,tol)
  x = x0;
  r = b - A * x0;
  p0 = r;
  s0 = A * p0;
  p1 = p0;
  s1 = s0;
  
  // .... first iteration 
  
    FloatVector p2 = p1;
    FloatVector s2 = s1;
    p1 = p0;
    s1 = s0;
    
    Float alpha = ( r * s1 ) / ( s1 * s1);
    x += alpha * p1;
    r -= alpha * s1;
    
    if ( r * r < tol*tol )
        break;
    
    p0 = s1;
    s0 = A * s1;
    
    Float beta1 = ( s0 * s1 ) / ( s1 * s1 );
    p0 -= beta1 * p1;
    s0 -= beta1 * s1;
            
  // ..... end first iteration 
  
  for iter=[1:maxit]
    
    FloatVector p2 = p1;
    FloatVector s2 = s1;
    p1 = p0;
    s1 = s0;
    
    Float alpha = ( r * s1 ) / ( s1 * s1);
    x += alpha * p1;
    r -= alpha * s1;
    
    if ( r * r < tol*tol )
        break;
    
    p0 = s1;
    s0 = A * s1;
    
    Float beta1 = ( s0 * s1 ) / ( s1 * s1 );
    p0 -= beta1 * p1;
    s0 -= beta1 * s1;
    
    if ( iter > 1 ) {
    
      Float beta2 = ( s0 * s2 ) / (s2 * s2 );
      p0 -= beta2 * p2;
      s0 -= beta2 * s2;
      
    }
    
  end
end*/  
































ResidualDescentMethod::ResidualDescentMethod( const LinearOperator& op )
: IterativeSolver(), A( op )
{
    this->max_iteration_count = op.getdimin();
    ResidualDescentMethod::check();
}

ResidualDescentMethod::~ResidualDescentMethod()
{
    ResidualDescentMethod::check();
}

void ResidualDescentMethod::check() const
{
    IterativeSolver::check();
    assert( A.getdimin() == A.getdimout() );
}

std::string ResidualDescentMethod::text() const
{
    return "Solver: Residual Descent Method";
}



// NOTE
// This method solves a possibly non-symmetric least squares system 
// by shifting the solution vector x within the span of the residual.


void ResidualDescentMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
    const int dimension = A.getdimin();

    /* Build up data */
    
    Float r_r = notanumber;

    FloatVector  r( dimension, 0. );
    
    // avoid repeated allocation of this temporary vector 
    FloatVector Ar( dimension, 0. );
    
    recent_iteration_count = 0;
    
    if( verbosity >= VerbosityLevel::verbose ) LOG << "START Residual iteration" << nl;
                
    while( recent_iteration_count < max_iteration_count )
    {
        
        /* Start / Restart CRM process */
        if( recent_iteration_count % x.getdimension() == 0 ) {
        
            r = b - A * x;

            r_r = r * r;
            
            assert( r_r >= 0. );
            
            if( verbosity >= VerbosityLevel::verbose )
                LOGPRINTF( "RESTARTED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(r_r), (long double)threshold );

        }

        bool continue_condition = r_r > threshold;
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( verbosity >= VerbosityLevel::verbose and print_condition )
            LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(r_r), (long double)threshold );
        
        /* If exit condition met, exit */
        
        if( not continue_condition ) 
            break;
            
        /* Perform iteration step */
        {

            Ar = A * r;
        
            Float  Ar_r = Ar * r;
            Float Ar_Ar  = Ar * Ar;
            assert( Ar_r >= 0. && Ar_Ar >= 0. );
            Float alpha = Ar_r / Ar_Ar;

            x += alpha * r;

            r -= alpha * Ar;
            
            r_r = r*r;

        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
    
    /* HOW DID WE FINISH ? */
    recent_deviation = r_r;
    
    if( verbosity >= VerbosityLevel::resultonly and print_modulo >= 0 ) 
        LOGPRINTF( "FINISHED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(recent_deviation), (long double)threshold );

    
}
  

  
  
  
/* 
 * 
 * Input: A, b, x
 *****************
 *
 * r := b - A x
 * Ar = A r
 * 
 * loop until exit condition
 * 
 *   rAr = r . Ar 
 *   rAtAr = Ar . Ar 
 *   
 *   alpha = rAr / rAtAr
 *   
 *   x += alpha r
 *   r -= alpha Ar
 *   
 *   
 */  

 































HerzogSoodhalterMethod::HerzogSoodhalterMethod( const LinearOperator& op )
: IterativeSolver(), A( op )
{
    this->max_iteration_count = op.getdimin();
    HerzogSoodhalterMethod::check();
}

HerzogSoodhalterMethod::~HerzogSoodhalterMethod()
{
    HerzogSoodhalterMethod::check();
}

void HerzogSoodhalterMethod::check() const
{
    IterativeSolver::check();
    assert( A.getdimin() == A.getdimout() );
}

std::string HerzogSoodhalterMethod::text() const
{
    return "Solver: Herzog-Soodhalter Method";
}



  
void HerzogSoodhalterMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
    const int dimension = A.getdimin();

    FloatVector v0( dimension, 0. );
    FloatVector v1( dimension, 0. );
    FloatVector w0( dimension, 0. );
    FloatVector w1( dimension, 0. );
    
    FloatVector vn( dimension, 0. );
    FloatVector wn( dimension, 0. );
    FloatVector  p( dimension, 0. );
    
    Float gamma = notanumber;
    Float eta = notanumber;
    
    Float s0 = notanumber;
    Float s1 = notanumber;
    Float c0 = notanumber;
    Float c1 = notanumber;
    
    recent_iteration_count = 0;

    while( recent_iteration_count < max_iteration_count ){
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );;
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and ( absolute(eta) < threshold );
        
        if( restart_condition or ( cpp_restart_before_finish and residual_seems_small ) ) {
            
            v0.zero();
            w0.zero();
            v1 = b - A * x;
            w1.zero();
            
            gamma = v1.norm();
            v1 /= gamma;
            
            assert( gamma > 0. );
            
            s0 = s1 = 0;
            c0 = c1 = 1;
            
            eta = gamma;
            
            if( verbosity >= VerbosityLevel::verbose ) {
                LOGPRINTF( "RESTARTED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) eta, (long double)threshold );
                LOGPRINTF( "NOTE Gamma: %.9Le Res: %.9Le\n", (long double)gamma, (long double)(b - A * x).norm() );
            }

        }
        
        bool residual_is_small = ( absolute(eta) < threshold );
        
        if( residual_is_small )
            break;

            
        {
            

            /*FloatVector*/ p = A * v1;
 
            Float delta = p * v1;
 
            /*FloatVector*/ vn = p - delta * v1 - gamma * v0;
 
            Float gamma_n = vn.norm();
            vn /= gamma_n;
            assert( gamma_n > 0. );
 
            Float alpha_0 = c1 * delta - c0 * s1 * gamma;
            assert( alpha_0 * alpha_0 + gamma_n * gamma_n );
            Float alpha_1 = std::sqrt( alpha_0 * alpha_0 + gamma_n * gamma_n );
            Float alpha_2 = s1 * delta + c0 * c1 * gamma;
            Float alpha_3 = s0 * gamma;
 
            assert( alpha_1 > 0. );

            Float cn = alpha_0 / alpha_1;
            Float sn = gamma_n / alpha_1;
            
            /*FloatVector*/ wn = ( v1 - alpha_2 * w1 - alpha_3 * w0 ) / alpha_1;
            x = x + cn * eta * wn;
 
            eta = - sn * eta;
            
            
            v0 = v1;
            w0 = w1;
            v1 = vn;
            w1 = wn;
            
            gamma = gamma_n;
            
            c0 = c1;
            s0 = s1;
            c1 = cn;
            s1 = sn;

        }

        
        
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        recent_deviation = eta; //( b - A * x ).norm_sq();
        
        if( verbosity >= VerbosityLevel::verbose and print_condition ) {
            LOGPRINTF( "INTERIM (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) recent_deviation, (long double)threshold );
            LOGPRINTF( "INTERIM Gamma: %.9Le Res: %.9Le\n", (long double)gamma, (long double)(b - A * x).norm() );
        }
        
        recent_iteration_count++;
        
    }
    
    /* HOW DID WE FINISH ? */
    
    recent_deviation = ( b - A * x ).norm_sq();
    
    if( verbosity >= VerbosityLevel::resultonly and print_modulo >= 0 ) 
        LOGPRINTF( "FINISHED (%d/%d) Residual: %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) sqrt(recent_deviation), (long double)threshold );

    
}
  

  


/*
 * 
 * v0 = 0
 * w0 = 0
 * v1 = b - A x
 * w1 = 0
 * 
 * gamma = norm(v1)
 * v1 = v1/gamma
 * 
 * s0 = 0
 * s1 = 0
 * c0 = 1
 * c1 = 1
 * 
 * eta = gamma
 * 
 * FOR j = 1 ..... 
 * 
 *  Av = A * v
 *  
 *  delta = Av * v
 *  
 *  v_new = Av - delta * v1 - gamma_1 v0
 *  
 *  gamma_new = norm( v_new )
 *  v_new /= gamma
 *  
 *  alpha_0 = c1 delta - c0 s1 gamma
 *  alpha_1 = sqrt( alpha_0 * alpha_0 + gamma_new * gamma_new )
 *  alpha_2 = s1 delta - c0 c1 gamma 
 *  alpha_3 = s0 gamma 
 * 
 *  c_new = alpha0 / alpha1
 *  s_new = gamma_new / alpha1
 *  
 *  w_new = ( v1 - alpha2 w1 - alpha3 w0 ) / alpha1
 *  x = x - c_new eta w_new
 *  
 *  eta = - s_new eta 
 * 
 * 
 * 
*/

