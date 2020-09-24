
#include "iterativesolver.hpp"

#include "../operators/floatvector.hpp"














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

void ConjugateGradientMethod::print( std::ostream& os ) const
{
    os << "Print Conjugate Gradient Method." << std::endl;
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
    
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = std::sqrt( r * r ) < tolerance;
        
        /* Start / Restart CRM process */
        if( restart_condition or residual_seems_small ) {
        
            LOG << "Begin Conjugate Gradient iteration";
        
            r = b - A * x;
            d = r;
            
            LOG << "starting with"
                << " r-norm: "    << std::sqrt( r * r ) 
                << " tolerance: " << tolerance;

        }

        bool residual_is_small = std::sqrt( r * r ) < tolerance;
        
        /* Print information if it is time too */
        bool do_print = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 ) or residual_is_small;
        if( do_print and verbosity >= VerbosityLevel::verbose ) {
            LOG << " # "  << recent_iteration_count << "/" << max_iteration_count
                << " $ "  << std::sqrt( r * r ) << "/" << tolerance;
        }

        /* If exit condition met, exit */
        if( residual_is_small ) 
            break;
            
        /* Perform iteration step */
        {
            
            Ad = A * d;
            
            Float rr_old = r * r;           assert( rr_old >= 0 ); if( std::sqrt(rr_old) < tolerance ){ break; }
            Float Ad_d  = Ad * d;
            Float alpha = rr_old / Ad_d;
            
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
    recent_deviation = std::sqrt(r * r);
    if( recent_deviation > tolerance ) {
            LOG << "CGM process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
        } else { 
            LOG << "CGM process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
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
        bool do_print = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 ) or residual_is_small;
        if( do_print and verbosity >= VerbosityLevel::verbose ) {
            LOG << " # "  << recent_iteration_count << "/" << max_iteration_count
                << " $ "  << std::sqrt( r * r ) << "/" << tolerance
                << " $ "  << std::sqrt( rAr ) << "/" << tolerance;
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
//             if( std::sqrt(Ad_Ad) < tolerance ) {
//                 LOG << "premature termination" << nl;
//                 return;
//             }
            
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
        
        bool do_print = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 ) or residual_is_small;
        if( do_print and verbosity >= VerbosityLevel::verbose ) {
            LOG << " # "  << recent_iteration_count << "/" << max_iteration_count
                << " $ "  << std::sqrt( r * r ) << "/" << tolerance;
        }

        
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
        
        bool do_print = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 ) or residual_is_small;
        if( do_print and verbosity >= VerbosityLevel::verbose ) {
            LOG << " # "  << recent_iteration_count << "/" << max_iteration_count
                << " $ "  << std::sqrt( r * r ) << "/" << tolerance;
        }

        
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

void PreconditionedConjugateResidualMethod::print( std::ostream& os ) const
{
    os << "Print Preconditioned Conjugate Residual Methop." << std::endl;
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
        
    LOG << "Begin Preconditioned Conjugate Residual iteration";// << std::endl;
    
    recent_iteration_count = 0;
    
    while(true)
    {
        
        /* Start / Restart PCRM process */
        if( recent_iteration_count % x.getdimension() == 0 ) {
        
            LOG << "Begin Preconditioned Conjugate Residual iteration";// << std::endl;
        
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
                
            }

            LOG << "starting with"
                      << " r-MAMsqnorm=" << rMAMr
                      << " r-Msqnorm="   << Mr * Mr  
                      << " r-sqnorm="    << r * r
                      ;//<< std::endl;

            LOG << "tolerance: " << tolerance;// << std::endl;

        }

        bool continue_condition = recent_iteration_count < max_iteration_count && rMAMr > tolerance;
        
        /* Print information if it is time too */
        if( recent_iteration_count % print_modulo == 0 or not continue_condition ) {
            LOG 
                << "#" << recent_iteration_count << "/" << max_iteration_count
                << " r-MAMsqnorm=" << rMAMr
                << " r-Msqnorm="   << ( r * ( M * r ) ) 
                ;//<< std::endl;
        }

        /* If exit condition met, exit */
        if( not continue_condition ) 
            break;
            
        /* Perform iteration step */
        {
            
            MAMp  = M *  AMp;
            AMAMp = A * MAMp;

            Float alpha = rMAMr / ( AMp * MAMp );

            x += alpha * Mp;

            r  -= alpha *   AMp; // not necessary, see also below
            Mr -= alpha *  MAMp;
            AMr -= alpha * AMAMp;

            Float newrMAMr = Mr * AMr;

            Float beta = newrMAMr / rMAMr;

            p =   r + beta *   p; // as such, completely useless
            Mp =  Mr + beta *  Mp;
            AMp = AMr + beta * AMp;

            rMAMr = newrMAMr;

        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
        

    /* HOW DID WE FINISH ? */
    recent_deviation = rMAMr;
    if( rMAMr > tolerance ) {
        LOG << "PCRM process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
    } else { 
        LOG << "PCRM process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;

    }

    
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

void MinimumResidualMethod::print( std::ostream& os ) const
{
    os << "Print MinimumResidualMethod." << std::endl;
}



  
void MinimumResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
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
    
    while( recent_iteration_count < max_iteration_count )
    {
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( recent_iteration_count % x.getdimension() == 0 );
        
        bool residual_seems_small = std::sqrt( rr ) < tolerance;
        
        /* Start / Restart MinimumResidualMethod process */
        if( restart_condition or residual_seems_small ) {
        
            if( verbosity >= VerbosityLevel::verbose ) 
            LOG << "Begin Minimal Residual iteration";// << std::endl;
        
            {
                
                r = b - A * x;
                p0 = r;
                s0 = A * p0;
                p1 = p0;
                s1 = s0;

                // .... first iteration 
            
                Float s1_s1 = s1 * s1;
                
                if( issmall( s1_s1 ) )
                    return;
                
                Float alpha = ( r * s1 ) / ( s1_s1);
                x += alpha * p1;
                r -= alpha * s1;
                
                rr = r * r;
                
                p0 = s1;
                s0 = A * s1;
                
                Float beta1 = ( s0 * s1 ) / ( s1_s1 );
                p0 -= beta1 * p1;
                s0 -= beta1 * s1;

            }
        
            if( verbosity >= VerbosityLevel::verbose ) {
                LOG << " # "  << recent_iteration_count << "/" << max_iteration_count
                    << " $ "  << std::sqrt( r * r ) << "/" << tolerance;
            }

        }

        bool residual_is_small = std::sqrt( rr ) < tolerance; 
        
        /* Print information if it is time too */
        bool do_print = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        if( do_print and verbosity >= VerbosityLevel::verbose ) 
        if( do_print or residual_is_small ) {
            LOG << "#" << recent_iteration_count << "/" << max_iteration_count
                << " r-sqnorm="  << rr 
                ;//<< std::endl;
        }

        /* If exit condition met, exit */
        if( residual_is_small ) 
            break;
            
            
        /* Perform iteration step */
        {
            
            p2 = p1; p1 = p0;
            s2 = s1; s1 = s0;
            
            Float s1_s1 = s1 * s1;
            
//             if( issmall( s1_s1 ) )
//                 break;

            Float alpha = ( r * s1 ) / ( s1_s1 );
            x += alpha * p1;
            r -= alpha * s1;
            
            rr = r * r;
            
            p0 = s1;
            s0 = A * s1;
            
            Float beta1 = ( s0 * s1 ) / ( s1 * s1 );
            p0 -= beta1 * p1;
            s0 -= beta1 * s1;
            
            Float s2_s2 = s2 * s2;
            
//             if( issmall( s2_s2 ) )
//                 break;

            Float beta2 = ( s0 * s2 ) / ( s2_s2 );
            p0 -= beta2 * p2;
            s0 -= beta2 * s2;
            
        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
    
    /* HOW DID WE FINISH ? */
    recent_deviation = std::sqrt(rr);
    if( verbosity >= VerbosityLevel::resultonly ) {
        if( recent_deviation > tolerance ) {
            LOG << "Minimum Residual process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
        } else { 
            LOG << "Minimum Residual process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
        }
    }
    
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

void ResidualDescentMethod::print( std::ostream& os ) const
{
    os << "Print Residual Descent Method." << std::endl;
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
    
    Float r_sqnorm = notanumber;

    FloatVector  r( dimension, 0. );
    
    // avoid repeated allocation of this temporary vector 
    FloatVector Ar( dimension, 0. );
    
    recent_iteration_count = 0;
    
    while(true)
    {
        
        /* Start / Restart CRM process */
        if( recent_iteration_count % x.getdimension() == 0 ) {
        
            LOG << "Begin Residual iteration";// << std::endl;
                
            r = b - A * x;

            r_sqnorm = r * r;
            
            assert( r_sqnorm >= 0. );
            
            LOG << "tolerance: " << tolerance;// << std::endl;

        }

        bool continue_condition = recent_iteration_count < max_iteration_count && r_sqnorm > tolerance;
        
        /* Print information if it is time too */
        if( recent_iteration_count % print_modulo == 0 or not continue_condition ) {
            LOG 
                << "#" << recent_iteration_count << "/" << max_iteration_count
                << " r-sqnorm=" << r_sqnorm 
                ;//<< std::endl;
        }

        /* If exit condition met, exit */
        if( not continue_condition ) 
            break;
            
        /* Perform iteration step */
        {

            Ar = A * r;
        
            Float  r_Asqnorm = r * Ar;
            Float Ar_sqnorm  = Ar * Ar;
            assert( r_Asqnorm >= 0. && Ar_sqnorm >= 0. );
            Float alpha = r_Asqnorm / Ar_sqnorm;

            x += alpha * r;

            r -= alpha * Ar;
            
            r_sqnorm = r*r;

        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
    
    /* HOW DID WE FINISH ? */
    recent_deviation = r_sqnorm;
    if( recent_deviation > tolerance ) {
        LOG << "Residual descent process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
    } else { 
        LOG << "Residual descent process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;

    }

    
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

void HerzogSoodhalterMethod::print( std::ostream& os ) const
{
    os << "Print Herzog-Soodhalter Method." << std::endl;
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
        
        
        bool restart_condition = (recent_iteration_count == 0);
        
        bool residual_seems_small = (eta < tolerance);
        
        if( restart_condition or residual_seems_small ) {
            
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
            
        }
        
        bool residual_is_small = (eta < tolerance);
        
        if( residual_is_small )
            break;

            
        {
            

            FloatVector p = A * v1;
 
            Float delta = p * v1;
 
            FloatVector vn = p - delta * v1 - gamma * v0;
 
            Float gamma_n = vn.norm();
            vn /= gamma_n;
            assert( gamma_n > 0. );
 
            Float alpha_0 = c1 * delta - c0 * s1 * gamma;
            Float alpha_1 = std::sqrt( alpha_0 * alpha_0 + gamma_n * gamma_n );
            Float alpha_2 = s1 * delta + c0 * c1 * gamma;
            Float alpha_3 = s0 * gamma;
 
            assert( alpha_1 > 0. );

            Float cn = alpha_0 / alpha_1;
            Float sn = gamma_n / alpha_1;
            
            FloatVector wn = ( v1 - alpha_2 * w1 - alpha_3 * w0 ) / alpha_1;
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

        
        
        
        bool do_print = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        if( do_print and verbosity >= VerbosityLevel::verbose ) {
            LOG << " # "  << recent_iteration_count << "/" << max_iteration_count
                << " $ "  << ( b - A * x ).norm() << "/" << tolerance;
        }
        
        
        recent_iteration_count++;
        
    }
    
    /* HOW DID WE FINISH ? */
    
    Float resnorm = ( b - A * x ).norm();
    recent_deviation = resnorm;
    
    if( recent_deviation > tolerance ) {
        LOG << "MINRES process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
    } else { 
        LOG << "MINRES process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
    }

    
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

