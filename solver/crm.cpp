
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
    
    while(true)
    {
        
        /* Start / Restart CRM process */
        if( recent_iteration_count % x.getdimension() == 0 ) {
        
            std::cout << "Begin Conjugate Residual iteration" << std::endl;
        
            iterationStart( x, b, r, d, Ar, Ad, rAr );
        
            std::cout << "starting with"
                      << " r-Asqnorm="   << r * Ar  
                      << " r-sqnorm="    << r * r 
                      << std::endl;
            std::cout << "tolerance: " << tolerance << std::endl;

        }

        bool continue_condition = recent_iteration_count < max_iteration_count && r * r > tolerance;
        
        /* Print information if it is time too */
        if( recent_iteration_count % print_modulo == 0 or not continue_condition ) {
            std::cout 
                << "#" << recent_iteration_count << "/" << max_iteration_count
                << " r-Asqnorm=" << rAr 
                << " r-sqnorm="  << r * r 
                << std::endl;
        }

        /* If exit condition met, exit */
        if( not continue_condition ) 
            break;
            
        /* Perform iteration step */
        iterationStep( x, r, d, Ar, Ad, rAr, AAd );
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
    
    /* HOW DID WE FINISH ? */
    if( rAr > tolerance ) {
        std::cout << "CRM process has failed.\n";
    } else { 
        std::cout << "CRM process has succeeded.\n";
    }

    recent_deviation = rAr;
    
}
  
  
  
void ConjugateResidualMethod::iterationStart( 
    const FloatVector& x, const FloatVector& b, 
    FloatVector& r, FloatVector& d, FloatVector& Ar, FloatVector& Ad,
    Float& rAr
) const {
    
    /* r = b - A x */
    r = b - A * x;

    /* d = r */
    d = A * r; //d.copydatafrom( r );

    /* Ar = A r */
    Ar = A * r; // A.apply( Ar, (const FloatVector&) r );

    /* Ad = A d */
    Ad = A * d; // A.apply( Ad, (const FloatVector&) d );

    /* rho is r.A.r */
    rAr = r * Ar;
    assert( rAr >= 0. );
    
}




void ConjugateResidualMethod::iterationStep( 
    FloatVector& x, 
    FloatVector& r, FloatVector& d, FloatVector& Ar, FloatVector& Ad,
    Float& rAr,
    FloatVector& AAd
) const {

    assert( rAr >= 0. );
    
    AAd = A * Ad;

    /*  alpha = r.A.r / d.AAd */
    Float Ad_Ad = Ad * Ad;
    Float alpha = rAr / Ad_Ad;

    assert( rAr >= 0. ); assert( Ad_Ad >= 0. );
    if( Ad_Ad < tolerance ) return;

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

  

  
void ConjugateResidualMethod::solve_robust( FloatVector& x, const FloatVector& b ) const
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
    
    while(true)
    {
        
        /* Start / Restart CRM process */
        if( recent_iteration_count % x.getdimension() == 0 ) {
        
            std::cout << "Begin Conjugate Residual iteration" << std::endl;
        
            r = b - A * x;
            d = A * r;
            
            std::cout << "starting with"
                      << " r-sqnorm="    << r * r 
                      << std::endl;
            std::cout << "tolerance: " << tolerance << std::endl;

        }

        bool continue_condition = recent_iteration_count < max_iteration_count && r * r > tolerance && d * d > tolerance;
        
        /* Print information if it is time too */
        if( recent_iteration_count % print_modulo == 0 or not continue_condition ) {
            std::cout 
                << "#" << recent_iteration_count << "/" << max_iteration_count
                << " r-sqnorm="  << r * r 
                << std::endl;
        }

        /* If exit condition met, exit */
        if( not continue_condition ) 
            break;
            
        /* Perform iteration step */
        {
            Ad = A * d;
            Float Ad_Ad = Ad * Ad; assert( Ad_Ad >= 0 ); if( Ad_Ad < tolerance ) break;
            Float Ad_r  = Ad * r;
            Float alpha = Ad_r / Ad_Ad;
            x = x + alpha * d;
            r = r - alpha * Ad;
            Float Ad_Ar = Ad * ( A * r );
            Float beta = Ad_Ar / Ad_Ad;
            d = r - beta * d;
        }
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
    
    /* HOW DID WE FINISH ? */
    if( r * r > tolerance ) {
        std::cout << "CRM process has failed.\n";
    } else { 
        std::cout << "CRM process has succeeded.\n";
    }

    recent_deviation = r * r;
    
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
