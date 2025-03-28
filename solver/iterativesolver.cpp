
#include <cmath>

#include <string>
#include <utility>

#include "iterativesolver.hpp"

#include "../operators/floatvector.hpp"


static const bool cpp_restart_on_full_dimension = false;
static const bool cpp_restart_before_finish     = false;










std::string ConjugateGradientMethod::text() const
{
    return "Solver: Conjugate Gradient Method";
}

PreconditionedConjugateGradientMethod::PreconditionedConjugateGradientMethod( const LinearOperator& op, const LinearOperator& M )
: IterativeSolver( op ), M( M )
{
    PreconditionedConjugateGradientMethod::check();
}

void PreconditionedConjugateGradientMethod::check() const
{
    IterativeSolver::check();
    M.check();
    assert( M.getdimin() == M.getdimout() );
    assert( A.getdimin() == M.getdimout() );
}

std::string PreconditionedConjugateGradientMethod::text() const
{
    return "Solver: Preconditioned Conjugate Gradient Method.";
}



std::string ConjugateResidualMethod::text() const
{
    return "Solver: Conjugate Residual Method";
}

PreconditionedConjugateResidualMethod::PreconditionedConjugateResidualMethod( const LinearOperator& op, const LinearOperator& M )
: IterativeSolver( op ), M( M )
{
    PreconditionedConjugateResidualMethod::check();
}

void PreconditionedConjugateResidualMethod::check() const
{
    IterativeSolver::check();
    M.check();
    assert( M.getdimin() == M.getdimout() );
    assert( A.getdimin() == M.getdimout() );
}

std::string PreconditionedConjugateResidualMethod::text() const
{
    return "Solver: Preconditioned Conjugate Residual Method.";
}


std::string HerzogSoodhalterMethod::text() const
{
    return "Solver: Herzog-Soodhalter Method";
}

std::string ResidualDescentMethod::text() const
{
    return "Solver: Residual Descent Method";
}

std::string MinimumResidualMethod::text() const
{
    return "Solver: Minimum Residual Method";
}










  

  
void ConjugateGradientMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
    /* Determine the print flags */

    const bool do_print_begin     = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_interim   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_restart   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_breakdown = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_warning   = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_finish    = verbosity >= VerbosityLevel::startandfinish;

    /* Build up data */

    const Float tolerance = precision * b.norm();
    
    const int dimension = A.getdimin();

    FloatVector  r( dimension, 0. );
    FloatVector  d( dimension, 0. );
    
    // avoid repeated allocation of these temporary vectors 
    FloatVector Ad( dimension, 0. );
    
    recent_iteration_count = 0;

    // Float sigma_min_sq = b * ( A * b ) / (b*b);
    
    if( do_print_begin ) 
        LOGPRINTF("(%d/%d)     BEGIN: %s\n", recent_iteration_count, max_iteration_count, "Conjugate Gradient iteration" );
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and absolute( r * r ) < tolerance*tolerance;
        
        /* Start / Restart CRM process */
        if( restart_condition or ( residual_seems_small and cpp_restart_before_finish ) ) UNLIKELY {
        
            r = b - A * x;
            d = r;

            if( do_print_restart )
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(r*r), (double)(safedouble)tolerance );
            
        }

        bool residual_is_small = absolute( r * r ) < tolerance*tolerance;
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( do_print_interim and print_condition )
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(r*r), (double)(safedouble)tolerance );
        
        /* If exit condition met, exit */
        
        if( residual_is_small ) 
            break;
            
        /* Perform iteration step */
        {
            
            Ad = A * d;
            
            Float r_r_old = r * r;
            Float Ad_d  = Ad * d;
            Float alpha = r_r_old / Ad_d;

            assert( Ad_d >= 0 );

            // sigma_min_sq = minimum( sigma_min_sq, Ad_d / (d*d) );
            // LOG << "@" << recent_iteration_count << " : " << std::sqrt( r_r_old / ( sigma_min_sq * (x*x) ) ) << " with eigenvalue bound " << std::sqrt(sigma_min_sq) << nl;

            bool numerator_is_unreasonable   = not std::isfinite(r_r_old) or r_r_old < 0.;
            bool denominator_is_unreasonable = not std::isfinite(Ad_d) or Ad_d < 0.;
            bool denominator_is_small        = std::sqrt(absolute(Ad_d)) < machine_epsilon;
            
            if( numerator_is_unreasonable or denominator_is_unreasonable ) {
                if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient energy is unreasonable with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)Ad_d );
                break;
            }
            
            if( denominator_is_small ) {
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(r_r_old), (double)(safedouble)tolerance );
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient energy is small with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)Ad_d );
                break;
            }
            
            x = x + alpha * d;
            r = r - alpha * Ad;
            
            Float r_r_new = r * r;
            Float beta = r_r_new / r_r_old;
            d = r + beta * d;
        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
    
    /* HOW DID WE FINISH ? */
    recent_deviation = absolute(r * r);
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, recent_deviation < tolerance ? "SUCCESS" : "FAIL", (double)(safedouble) std::sqrt(recent_deviation), (double)(safedouble)tolerance ); 
        
    
}



void PreconditionedConjugateGradientMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );

    /* Determine the print flags */

    const bool do_print_begin     = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_interim   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_restart   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_breakdown = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_warning   = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_finish    = verbosity >= VerbosityLevel::startandfinish;

    /* Build up data */
    
    const Float tolerance = precision * b.norm();
    
    const int dimension = A.getdimin();

    Float rMr = notanumber;
    
    FloatVector   r( dimension, 0. );
    FloatVector   p( dimension, 0. ); 
    
    // avoid repeated allocation of these temporary vectors 
    FloatVector  Mr( dimension, 0. );
    FloatVector  Ap( dimension, 0. );    
    
    recent_iteration_count = 0;
    
    if( do_print_begin ) 
        LOGPRINTF("(%d/%d)     BEGIN: %s\n", recent_iteration_count, max_iteration_count, "Preconditioned Conjugate Gradient iteration" );
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and absolute( rMr ) < tolerance*tolerance;
        
        /* Start / Restart PCRM process */
        if( restart_condition or ( residual_seems_small and cpp_restart_before_finish ) ) UNLIKELY {
        
            {
                
                /* x = initial guess */
                
                /* r = b - A x */
                /* p = M * r */
                
                r = b - A * x;
                Mr = M * r;
                p = Mr;

                /* compute r * Mr */
                rMr = r * Mr;

                if( do_print_restart )
                    LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(rMr), (double)(safedouble)tolerance );

            }

        }

        /* If exit condition met, exit */

        bool residual_is_small = absolute( rMr ) < tolerance*tolerance;
        
        if( residual_is_small ) 
            break;
            
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );

        if( do_print_interim and print_condition ) 
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(rMr), (double)(safedouble)tolerance );

        /* Perform iteration step */
        {
            
            Ap = A * p;
            
            Float p_Ap = p * Ap;
            Float alpha = rMr / p_Ap;

            bool numerator_is_unreasonable   = not std::isfinite(p_Ap) or p_Ap < 0.;
            bool denominator_is_unreasonable = not std::isfinite(p_Ap) or p_Ap < 0.;
            bool denominator_is_small        = std::sqrt(absolute(p_Ap)) < machine_epsilon;
            
            if( numerator_is_unreasonable or denominator_is_unreasonable ) {
                if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient energy is unreasonable with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)p_Ap );
                break;
            }
            
            if( denominator_is_small ) {
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(r * Mr), (double)(safedouble)tolerance );
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient energy is small with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)p_Ap );
                break;
            }
            
            x += alpha * p;
            r -= alpha * Ap;
            
            Mr = M * r;

            Float new_rMr = r * Mr;

            Float beta = new_rMr / rMr;

            p =   Mr + beta * p; 

            rMr = new_rMr;

        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
        

    /* HOW DID WE FINISH ? */

    recent_deviation = rMr;
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, recent_deviation < tolerance ? "SUCCESS" : "FAIL", (double)(safedouble) std::sqrt(recent_deviation), (double)(safedouble)tolerance );

    
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
    
    /* Determine the print flags */

    const bool do_print_begin     = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_interim   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_restart   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_breakdown = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_warning   = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_finish    = verbosity >= VerbosityLevel::startandfinish;

    /* Build up data */
    
    const Float tolerance = precision * b.norm();
    
    const int dimension = A.getdimin();

    Float rAr = notanumber;

    FloatVector  r( dimension, 0. );
    FloatVector  d( dimension, 0. );
    
    FloatVector Ad( dimension, 0. );
    FloatVector Ar( dimension, 0. );
    
    FloatVector AAd( dimension, 0. );
    
    recent_iteration_count = 0;
    
    if( do_print_begin ) 
        LOGPRINTF("(%d/%d)     BEGIN: %s\n", recent_iteration_count, max_iteration_count, "Conjugate Residual (explicit) iteration" );
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        /* Start / Restart CRM process */
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and ( absolute( r * r ) < tolerance*tolerance or absolute( rAr ) < tolerance*tolerance );
        
        if( restart_condition or ( residual_seems_small and cpp_restart_before_finish ) ) UNLIKELY {
        
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

            if( do_print_restart )
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(rAr), (double)(safedouble)tolerance );
            

            
        }

        bool residual_is_small = absolute( r * r ) < tolerance*tolerance or absolute( rAr ) < tolerance*tolerance; 
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( do_print_interim and print_condition )
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(rAr), (double)(safedouble)tolerance ); 
        
        /* If exit condition met, exit */
        if( residual_is_small ) 
            break;
            
        
        /* Perform iteration step */
        {

            if( rAr < 0. ) 
            // LOG << rAr << nl;
            assert( rAr >= 0. );
            
            AAd = A * Ad;

            /*  alpha = r.A.r / d.AAd */
            Float Ad_Ad = Ad * Ad;
            Float alpha = rAr / Ad_Ad;

            bool numerator_is_unreasonable   = not std::isfinite(rAr) or rAr < 0.;
            bool denominator_is_unreasonable = not std::isfinite(Ad_Ad) or Ad_Ad < 0.;
            bool denominator_is_small        = std::sqrt(absolute(Ad_Ad)) < machine_epsilon;
            
            if( numerator_is_unreasonable or denominator_is_unreasonable ) {
                if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient double energy is unreasonable with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)Ad_Ad );
                break;
            }
            
            if( denominator_is_small ) {
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(rAr), (double)(safedouble)tolerance ); 
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient double energy is small with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)Ad_Ad );
                break;
            }
            
                        
            x += alpha * d;

            r -= alpha * Ad; //= b - A * x;

            Ar -= alpha * AAd;

            Float tau = rAr;
            rAr = r * Ar; //r * ( A * r ); //r * Ar;
            Float beta = rAr / tau;
            
            if( rAr < 0. ) 
            {
                if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Residual energy is unreasonable with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)rAr );
                // rAr = 0.; 
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
    recent_deviation = ( b - A * x ).norm_sq();
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, recent_deviation < tolerance ? "SUCCESS" : "FAIL", (double)(safedouble) std::sqrt(recent_deviation), (double)(safedouble)tolerance ); 

    
}

void ConjugateResidualMethod::solve_robust( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
    /* Determine the print flags */

    const bool do_print_begin     = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_interim   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_restart   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_breakdown = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_warning   = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_finish    = verbosity >= VerbosityLevel::startandfinish;

    /* Build up data */
    
    const Float tolerance = precision * b.norm();
    
    const int dimension = A.getdimin();

    FloatVector  r( dimension, 0. );
    FloatVector  d( dimension, 0. );
    FloatVector Ar( dimension, 0. );
    FloatVector Ad( dimension, 0. );
    FloatVector  p( dimension, 0. );
    
    recent_iteration_count = 0;
    
    if( do_print_begin ) 
        LOGPRINTF("(%d/%d)     BEGIN: %s\n", recent_iteration_count, max_iteration_count, "Conjugate Residual (robust) iteration" );
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and ( absolute( r * r ) < tolerance*tolerance or absolute( r * Ar ) < tolerance*tolerance );
        // first criterion is not in fast 

        if( restart_condition or ( residual_seems_small and cpp_restart_before_finish ) ) UNLIKELY {
        
            r  = b - A * x;
            d  = r;
            
            Ar = A * r;
            Ad = Ar;

            // fast and explicit: Ar_r = Ar * r;

            if( do_print_restart )
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt( r * Ar ), (double)(safedouble)tolerance );

        }

        bool residual_is_small = absolute( r * r ) < tolerance*tolerance or absolute( r * Ar ) < tolerance*tolerance; 
        // first criterion is in explicit but not in fast (for speed reasons)
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( do_print_interim and print_condition )
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt( Ar * r ), (double)(safedouble)tolerance ); 
        
        /* If exit condition met, exit */
        
        if( residual_is_small ) 
            break;
            
        {
            p = A * Ad;
            Float Ad_Ad = Ad * Ad; Assert( Ad_Ad >= 0, A, Ad_Ad, Ad, A * Ad ); 
            
            Float Ad_r  = Ad * r; // fast/explicit uses and maintains Ar_r
            Float alpha = Ad_r / Ad_Ad;
            
            bool numerator_is_unreasonable   = not std::isfinite(Ad_r) or Ad_r < 0.;
            bool denominator_is_unreasonable = not std::isfinite(Ad_Ad) or Ad_Ad < 0.;
            bool numerator_is_small          = std::sqrt(absolute(Ad_r)) < machine_epsilon;
            bool denominator_is_small        = std::sqrt(absolute(Ad_Ad)) < machine_epsilon;
            
            if( numerator_is_unreasonable or denominator_is_unreasonable ) {
                if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient double energy is unreasonable with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)Ad_Ad );
                break;
            }
            
            if( numerator_is_small or denominator_is_small ) {
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt( Ar * r ), (double)(safedouble)tolerance ); 
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient double energy is small with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)Ad_Ad );
                break;
            }
            
            x  =  x + alpha * d;
            r  =  r - alpha * Ad;
            Ar = Ar - alpha * p;
            
            Float Ar_r_new = Ar * r;
            Float beta = Ar_r_new / Ad_r; // fast uses and maintains Ar_r 
            
            // assert( Ar_r_new >= 0. );

            d  =  r + beta *  d;
            Ad = Ar + beta * Ad;
            
        }
        
        recent_iteration_count++;
    }
    
    
    recent_deviation = absolute( r * Ar );
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, recent_deviation < tolerance ? "SUCCESS" : "FAIL", (double)(safedouble) std::sqrt(recent_deviation), (double)(safedouble)tolerance ); 
    
}

void ConjugateResidualMethod::solve_fast( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
    /* Determine the print flags */

    const bool do_print_begin     = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_interim   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_restart   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_breakdown = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_warning   = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_finish    = verbosity >= VerbosityLevel::startandfinish;

    /* Build up data */

    const Float tolerance = precision * b.norm();
    
    const int dimension = A.getdimin();

    FloatVector  r( dimension, 0. );
    FloatVector  d( dimension, 0. );
    FloatVector Ar( dimension, 0. );
    FloatVector Ad( dimension, 0. );
    FloatVector  p( dimension, 0. );
    
    recent_iteration_count = 0;
    
    Float Ar_r = notanumber;
    
    if( do_print_begin ) 
        LOGPRINTF("(%d/%d)     BEGIN: %s\n", recent_iteration_count, max_iteration_count, "Conjugate Residual (fast) iteration" );
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and absolute( Ar_r ) < tolerance*tolerance;
        
        if( restart_condition or ( residual_seems_small and cpp_restart_before_finish ) ) UNLIKELY {
        
            r  = b - A * x;
            d  = r;
            
            Ar = A * r;
            Ad = Ar;
            
            Ar_r = Ar * r;

            if( do_print_restart )
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(Ar_r), (double)(safedouble)tolerance );

        }

        bool residual_is_small = absolute( Ar * r ) < tolerance*tolerance; 
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( do_print_interim and print_condition )
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt( Ar * r ), (double)(safedouble)tolerance );
        
        /* If exit condition met, exit */
        
        if( residual_is_small ) 
            break;
            
        {
            
            p = A * Ad;
            Float Ad_Ad = Ad * Ad; assert( Ad_Ad >= 0 ); 
            
            Float alpha = Ar_r / Ad_Ad;
            
            bool numerator_is_unreasonable   = not std::isfinite(Ar_r) or Ar_r < 0.;
            bool denominator_is_unreasonable = not std::isfinite(Ad_Ad) or Ad_Ad < 0.;
            bool numerator_is_small          = std::sqrt(absolute(Ar_r)) < machine_epsilon;
            bool denominator_is_small        = std::sqrt(absolute(Ad_Ad)) < machine_epsilon;
            
            if( numerator_is_unreasonable or denominator_is_unreasonable ) {
                if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient double energy is unreasonable with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)Ad_Ad );
                break;
            }
            
            if( numerator_is_small or denominator_is_small ) {
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt( Ar * r ), (double)(safedouble)tolerance );
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient double energy is small with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)Ad_Ad );
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
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, recent_deviation < tolerance ? "SUCCESS" : "FAIL", (double)(safedouble) std::sqrt(recent_deviation), (double)(safedouble)tolerance );
    
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

 
 


























void PreconditionedConjugateResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );

    /* Determine the print flags */

    const bool do_print_begin     = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_interim   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_restart   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_breakdown = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_warning   = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_finish    = verbosity >= VerbosityLevel::startandfinish;

    /* Build up data */
    
    const Float tolerance = precision * b.norm();
    
    const int dimension = A.getdimin();

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
    
    if( do_print_begin ) 
        LOGPRINTF("(%d/%d)     BEGIN: %s\n", recent_iteration_count, max_iteration_count, "Preconditioned Conjugate Residual iteration" );
        
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and absolute( rMAMr ) < tolerance*tolerance;
        
        /* Start / Restart PCRM process */
        if( restart_condition or ( residual_seems_small and cpp_restart_before_finish ) ) UNLIKELY {
        
            {
                
                r = b - A * x;
                p = M * r;
                
                Mr = M * r;
                Mp = M * p;

                AMr = A * Mr;
                AMp = A * Mp;

                /* rho is Mr.A.Mr */
                rMAMr = Mr * AMr;
                
                if( do_print_restart )
                    LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(rMAMr), (double)(safedouble)tolerance );

            }

        }

        /* If exit condition met, exit */

        bool residual_is_small = absolute( rMAMr ) < tolerance*tolerance;
        
        if( residual_is_small ) 
            break;
            
        /* Print information */
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );;

        if( do_print_interim and print_condition ) 
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(rMAMr), (double)(safedouble)tolerance );

        /* Perform iteration step */
        {
            
            MAMp  = M *  AMp;
            AMAMp = A * MAMp;

            Float AMp_MAMp = AMp * MAMp;
            Float alpha = rMAMr / ( AMp_MAMp );

            bool numerator_is_unreasonable   = not std::isfinite(rMAMr) or rMAMr < 0.;
            bool denominator_is_unreasonable = not std::isfinite(AMp_MAMp) or AMp_MAMp < 0.;
            bool denominator_is_small        = std::sqrt(absolute(AMp_MAMp)) < machine_epsilon;
            
            if( numerator_is_unreasonable or denominator_is_unreasonable ) {
                if( do_print_breakdown ) LOGPRINTF( "(%d/%d) BREAKDOWN: Gradient double energy is unreasonable with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)AMp_MAMp );
                break;
            }
            
            if( denominator_is_small ) {
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(Mr * AMr), (double)(safedouble)tolerance );
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Gradient double energy is small with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)AMp_MAMp );
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
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, recent_deviation < tolerance ? "SUCCESS" : "FAIL", (double)(safedouble) std::sqrt(recent_deviation), (double)(safedouble)tolerance );

    
}
  







































  
void MinimumResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    assert( x.is_finite() and b.is_finite() );
    
    /* Determine the print flags */

    const bool do_print_begin     = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_interim   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_restart   = verbosity >= VerbosityLevel::verbose;
    // const bool do_print_breakdown = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_warning   = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_finish    = verbosity >= VerbosityLevel::startandfinish;

    /* Build up data */
    
    const Float tolerance = precision * b.norm();
    
    const int dimension = A.getdimin();

    Float r_r = notanumber;
    
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
    
    if( do_print_begin ) 
        LOGPRINTF("(%d/%d)     BEGIN: %s\n", recent_iteration_count, max_iteration_count, "Minimal Residual iteration" );

    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and absolute( r_r ) < tolerance*tolerance;
        
        /* Start / Restart MinimumResidualMethod process */
        if( restart_condition or ( residual_seems_small and cpp_restart_before_finish ) ) UNLIKELY {
        
            r = b - A * x;
            
            r_r = r * r;
            
            assert( std::isfinite(r_r) );
            if( r_r < tolerance*tolerance )
                break;
            
            p0 = r;
            s0 = A * r;
            
            
            Float s0_s0 = s0 * s0;
            
            assert( std::isfinite(s0_s0) );
            if( s0_s0 < tolerance*tolerance )
                break;
            
            p0 /= std::sqrt(s0_s0);
            s0 /= std::sqrt(s0_s0);
            
            Float alpha0 = r * s0;
            
            x = x + alpha0 * p0;
            r = r - alpha0 * s0;
            
            r_r = r * r;
            
            assert( std::isfinite(r_r) );
            if( r_r < tolerance*tolerance )
                break;
            
            
            p1 = r; 
            s1 = A * r;
            
            Float beta0 = s0 * s1;
            
            p1 = p1 - beta0 * p0;
            s1 = s1 - beta0 * s0;
            
            Float s1_s1 = s1 * s1;
            
            assert( std::isfinite(s1_s1) );
            if( s1_s1 < tolerance*tolerance )
                break;
            
            p1 /= std::sqrt(s1_s1);
            s1 /= std::sqrt(s1_s1);
            
            Float alpha1 = r * s1;
            
            x = x + alpha1 * p1;
            r = r - alpha1 * s1;
            
            r_r = r * r;
            
//             recent_alpha = alpha1;

            if( do_print_restart )
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(r_r), (double)(safedouble)tolerance );

            
        }

        bool residual_is_small = absolute( r_r ) < tolerance*tolerance; 
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( do_print_interim and print_condition )
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(r_r), (double)(safedouble)tolerance );
        
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
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(r_r), (double)(safedouble) tolerance );
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Norm of search direction is small with %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(s2_s2) );
                if( do_print_warning ) LOGPRINTF( "(%d/%d)   WARNING: Machine Epsilon: %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) machine_epsilon );
                // break;
            }
            
            Float mysqrt = std::sqrt(s2_s2);
            p2 /= mysqrt;
            s2 /= mysqrt;
            
            Float alpha2 = r * s2;
            
            x = x + alpha2 * p2;
            r = r - alpha2 * s2;
            
//             r = b - A * x;
            
            r_r = r * r;
            
//             recent_alpha = alpha2;
            
            assert( std::isfinite(r_r) );
            if( std::sqrt(r_r) < tolerance )
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
    // recent_deviation = absolute(r_r);
    recent_deviation = absolute( ( b - A * x ).norm_sq() );
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, recent_deviation < tolerance ? "SUCCESS" : "FAIL", (double)(safedouble) std::sqrt(recent_deviation), (double)(safedouble)tolerance );
    
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
    
    if( r * r < tol*tol )
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
    
    if( r * r < tol*tol )
        break;
    
    p0 = s1;
    s0 = A * s1;
    
    Float beta1 = ( s0 * s1 ) / ( s1 * s1 );
    p0 -= beta1 * p1;
    s0 -= beta1 * s1;
    
    if( iter > 1 ) {
    
      Float beta2 = ( s0 * s2 ) / (s2 * s2 );
      p0 -= beta2 * p2;
      s0 -= beta2 * s2;
      
    }
    
  end
end*/  































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
    
    /* Determine the print flags */

    const bool do_print_begin     = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_interim   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_restart   = verbosity >= VerbosityLevel::verbose;
    // const bool do_print_breakdown = verbosity >= VerbosityLevel::startandfinish;
    // const bool do_print_warning   = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_finish    = verbosity >= VerbosityLevel::startandfinish;

    /* Build up data */
    
    const Float tolerance = precision * b.norm();
    
    const int dimension = A.getdimin();

    Float r_r = notanumber;

    FloatVector  r( dimension, 0. );
    
    // avoid repeated allocation of this temporary vector 
    FloatVector Ar( dimension, 0. );
    
    recent_iteration_count = 0;
    
    if( do_print_begin ) 
        LOGPRINTF("(%d/%d)     BEGIN: %s\n", recent_iteration_count, max_iteration_count, "Residual Descent iteration" );
                
    while( recent_iteration_count < max_iteration_count )
    {
        
        /* Start / Restart CRM process */
        if( recent_iteration_count % x.getdimension() == 0 ) {
        
            r = b - A * x;

            r_r = r * r;
            
            assert( r_r >= 0. );
            
            if( do_print_restart )
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(r_r), (double)(safedouble)tolerance );

        }

        bool continue_condition = r_r > tolerance;
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( do_print_interim and print_condition )
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) std::sqrt(r_r), (double)(safedouble)tolerance );
        
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
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, recent_deviation < tolerance ? "SUCCESS" : "FAIL", (double)(safedouble) std::sqrt(recent_deviation), (double)(safedouble)tolerance );

    
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

 





























  
void HerzogSoodhalterMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
    /* Determine the print flags */

    const bool do_print_begin     = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_interim   = verbosity >= VerbosityLevel::verbose;
    const bool do_print_restart   = verbosity >= VerbosityLevel::verbose;
    // const bool do_print_breakdown = verbosity >= VerbosityLevel::startandfinish;
    // const bool do_print_warning   = verbosity >= VerbosityLevel::startandfinish;
    const bool do_print_finish    = verbosity >= VerbosityLevel::startandfinish;

    /* Build up data */
    
    const Float tolerance = precision * b.norm();
    
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

    if( do_print_begin ) 
        LOGPRINTF("(%d/%d)     BEGIN: %s\n", recent_iteration_count, max_iteration_count, "Herzog-Soodhalter iteration" );

                
    while( recent_iteration_count < max_iteration_count ){
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cpp_restart_on_full_dimension and recent_iteration_count % x.getdimension() == 0 );;
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and ( absolute(eta) < tolerance );
        
        if( restart_condition or ( residual_seems_small and cpp_restart_before_finish ) ) UNLIKELY {
            
            // 1 -- 2 
            v0.zero();
            w0.zero();
            v1 = b - A * x;
            w1.zero();
            
            // 3
            gamma = v1.norm();
            v1 /= gamma;
            
            assert( gamma > 0. );
            
            // 7
            s0 = s1 = 0.;
            c0 = c1 = 1.;
            
            eta = gamma;
            
            if( do_print_restart ) {
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) eta, (double)(safedouble)tolerance );
                LOGPRINTF( "(%d/%d)            Gamma: %.9le Res: %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)gamma, (double)(safedouble)(b - A * x).norm() );
            }

        }
        
        bool residual_is_small = ( absolute(eta) < tolerance );
        
        if( residual_is_small )
            break;

            
        {
            

            // 9 
            /*FloatVector*/ p = A * v1;
 
            Float delta = p * v1;
 
            // 10 
            /*FloatVector*/ vn = p - delta * v1 - gamma * v0;
 
            // 12 -- 13
            Float gamma_n = vn.norm();
            vn /= gamma_n;
            assert( gamma_n > 0. );
 
            // 15 - 18
            Float alpha_0 = c1 * delta - c0 * s1 * gamma;
            assert( alpha_0 * alpha_0 + gamma_n * gamma_n > 0. );
            Float alpha_1 = std::sqrt( alpha_0 * alpha_0 + gamma_n * gamma_n );
            Float alpha_2 = s1 * delta + c0 * c1 * gamma;
            Float alpha_3 = s0 * gamma;
 
            assert( alpha_1 > 0. );

            // 19 
            Float cn = alpha_0 / alpha_1;
            Float sn = gamma_n / alpha_1;
            
            // 23 -- 24
            /*FloatVector*/ wn = ( v1 - alpha_2 * w1 - alpha_3 * w0 ) / alpha_1;
            x = x + cn * eta * wn;
 
            // 27 
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

        recent_deviation = eta;
        
        
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( do_print_interim and print_condition ) {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) recent_deviation, (double)(safedouble)tolerance );
            LOGPRINTF( "(%d/%d)   INTERIM: Gamma: %.9le Eta: %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)gamma, (double)(safedouble)eta );
        }
        
        recent_iteration_count++;
        
    }
    
    /* HOW DID WE FINISH ? */
    
    recent_deviation = ( b - A * x ).norm_sq();
    
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, recent_deviation < tolerance ? "SUCCESS" : "FAIL", (double)(safedouble) std::sqrt(recent_deviation), (double)(safedouble)tolerance );

    
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
 *  alpha_1 = std::sqrt( alpha_0 * alpha_0 + gamma_new * gamma_new )
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

