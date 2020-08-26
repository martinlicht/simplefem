
#include "crm.hpp"

#include "../operators/floatvector.hpp"



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

void ConjugateResidualMethod::print( std::ostream& os ) const
{
    os << "Print Conjugate Residual Method." << std::endl;
}



  
void ConjugateResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
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
    
    // avoid repeated allocation of these temporary vectors 
    FloatVector AAd( dimension, 0. );
    
    recent_iteration_count = 0;
    
    while( recent_iteration_count < max_iteration_count )
    {
        
        /* Start / Restart CRM process */
        bool restart_condition = ( recent_iteration_count == 0 ) or ( recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = std::sqrt( r * r ) < tolerance or std::sqrt( rAr ) < tolerance;
        
        if( restart_condition or residual_seems_small ) {
        
            LOG << "Begin Conjugate Residual iteration";// << std::endl;
        
            {
    
                /* r = b - A x */
                r = b - A * x;

                /* d = r */
                d = r; //d.copydatafrom( r );

                /* Ar = A r */
                Ar = A * r; // A.apply( Ar, (const FloatVector&) r );

                /* Ad = A d */
                Ad = A * d; // A.apply( Ad, (const FloatVector&) d );

                /* rho is r.A.r */
                rAr = r * Ar;
                assert( rAr >= 0. );
                
            }
        
            LOG << "starting with"
                << " r-Anorm: "   << std::sqrt( r * Ar )   
                << " r-snorm: "    << std::sqrt( r *  r ) 
                << " tolerance: " << tolerance;// << std::endl;

        }

        bool residual_is_small = std::sqrt( r * r ) < tolerance or std::sqrt( rAr ) < tolerance; 
        
        /* Print information if it is time too */
        if( recent_iteration_count % print_modulo == 0 or residual_is_small ) {
            LOG << "#" << recent_iteration_count << "/" << max_iteration_count
                << " r-Anorm: " << std::sqrt( rAr )   
                << " r-snorm: " << std::sqrt( r *  r );
        }

        /* If exit condition met, exit */
        if( residual_is_small ) 
            break;
            
        /* Perform iteration step */
        {

            assert( rAr >= 0. );
            
            AAd = A * Ad;

            /*  alpha = r.A.r / d.AAd */
            Float Ad_Ad = Ad * Ad;
            Float alpha = rAr / Ad_Ad;

            assert( rAr >= 0. ); assert( Ad_Ad >= 0. );
            if( std::sqrt(Ad_Ad) < tolerance ) {
                LOG << "premature termination" << nl;
                return;
            }
            
            x += alpha * d;

            r -= alpha * Ad; //= b - A * x;

            Ar -= alpha * AAd;

            Float tau = rAr;
            rAr = r * ( A * r ); //r * ( A * r ); //r * Ar;
            Float beta = rAr / tau;
            assert( rAr >= 0. );
            
            d = r + beta * d;

            Ad = Ar + beta * Ad;
            
        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
    
    /* HOW DID WE FINISH ? */
    recent_deviation = std::sqrt( rAr );
    if( std::sqrt( r * r ) > tolerance and std::sqrt( rAr ) > tolerance ) {
            LOG << "CRM process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
        } else { 
            LOG << "CRM process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
    }

    
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
    
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = std::sqrt( r * r ) < tolerance;
        
        if( restart_condition ) {
        
            if( verbosity >= VerbosityLevel::verbose ) LOG << "Begin Conjugate Residual iteration";// << std::endl;
        
            r  = b - A * x;
            d  = r;
            
            Ad = A * d;
            Ar = Ad;
            
//             if( verbosity >= VerbosityLevel::verbose ) 
//                 LOG << "starting with"
//                     << " r-sqnorm="    << r * r 
//                     ;//<< std::endl;
//             if( verbosity >= VerbosityLevel::verbose ) 
//                 LOG << "tolerance: " << tolerance;// << std::endl;

        }

        bool residual_is_small = std::sqrt( r * r ) < tolerance; 
        
//         /* Print information if it is time too */
//         if( verbosity >= VerbosityLevel::verbose ) 
//         if( recent_iteration_count % print_modulo == 0 or residual_is_small ) {
//             LOG 
//                 << "#" << recent_iteration_count << "/" << max_iteration_count
//                 << " r-sqnorm="  << r * r 
//                 ;//<< std::endl;
//         }

        
        if( residual_is_small ) 
            break;
            
        {
            p = A * Ad;
            Float Ad_Ad = Ad * Ad; assert( Ad_Ad >= 0 ); 
            
            Float Ad_r  = Ad * r;
            Float alpha = Ad_r / Ad_Ad;
            
            x  =  x + alpha * d;
            r  =  r - alpha * Ad;
            Ar = Ar - alpha * p;
            
            Float Ar_r_new = Ar * r;
            Float beta = Ar_r_new / Ad_r;
            
            d  =  r + beta *  d;
            Ad = Ar + beta * Ad;
            
        }
        
        recent_iteration_count++;
    }
    
    
    recent_deviation = std::sqrt( r * r );

    if( verbosity >= VerbosityLevel::resultonly ) {
        if( recent_deviation > tolerance ) {
            LOG << "CRM process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
        } else { 
            LOG << "CRM process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
        }
    }
    
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
    
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = std::sqrt( r * r ) < tolerance;
        
        if( restart_condition ) {
        
            if( verbosity >= VerbosityLevel::verbose ) LOG << "Begin Conjugate Residual iteration";// << std::endl;
        
            r  = b - A * x;
            d  = r;
            
            Ad = A * d;
            Ar = Ad;
            
            Ar_r = Ar * r;
//             if( verbosity >= VerbosityLevel::verbose ) 
//                 LOG << "starting with"
//                     << " r-sqnorm="    << r * r 
//                     ;//<< std::endl;
//             if( verbosity >= VerbosityLevel::verbose ) 
//                 LOG << "tolerance: " << tolerance;// << std::endl;

        }

        bool residual_is_small = std::sqrt( r * r ) < tolerance; 
        
//         /* Print information if it is time too */
//         if( verbosity >= VerbosityLevel::verbose ) 
//         if( recent_iteration_count % print_modulo == 0 or residual_is_small ) {
//             LOG 
//                 << "#" << recent_iteration_count << "/" << max_iteration_count
//                 << " r-sqnorm="  << r * r 
//                 ;//<< std::endl;
//         }

        
        if( residual_is_small ) 
            break;
            
        {
            
            p = A * Ad;
            Float Ad_Ad = Ad * Ad; assert( Ad_Ad >= 0 ); 
            
            Float alpha = Ar_r / Ad_Ad;
            
            x  =  x + alpha * d;
            r  =  r - alpha * Ad;
            Ar = Ar - alpha * p;
            
            Float Ar_r_new = Ar * r;
            Float beta = Ar_r_new / Ar_r;
            
            d  =  r + beta *  d;
            Ad = Ar + beta * Ad;
            
            Ar_r = Ar_r_new;
            
        }
        
        recent_iteration_count++;
    }
    
    
    recent_deviation = std::sqrt( r * r );

    if( verbosity >= VerbosityLevel::resultonly ) {
        if( recent_deviation > tolerance ) {
            LOG << "CRM process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
        } else { 
            LOG << "CRM process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
        }
    }
    
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

 
 
 
 
 
 
 
 
 
