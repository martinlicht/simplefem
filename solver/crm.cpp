
#include "crm.hpp"

#include "../operators/floatvector.hpp"



ConjugateResidualMethod::ConjugateResidualMethod( const LinearOperator& op )
: IterativeSolver(), A( op )
{
    this->max_iteration_count = op.getdimin();
    ConjugateResidualMethod::check();
}

ConjugateResidualMethod::~ConjugateResidualMethod()
{}

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
    
    Float r_Anorm = notanumber;

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
        
            iterationStart( x, b, r, d, Ar, Ad, r_Anorm );
        
            std::cout << "tolerance: " << tolerance << std::endl;

        }

        bool continue_condition = recent_iteration_count < max_iteration_count && r * r > tolerance;
        
        /* Print information if it is time too */
        if( recent_iteration_count % print_modulo == 0 or not continue_condition ) {
            std::cout 
                << "#" << recent_iteration_count << "/" << max_iteration_count
                << " r-Asqnorm=" << r_Anorm 
                << " r-sqnorm="  << r * r 
                << std::endl;
        }

        /* If exit condition met, exit */
        if( not continue_condition ) 
            break;
            
        /* Perform iteration step */
        iterationStep( x, r, d, Ar, Ad, r_Anorm, AAd );
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
    
    /* HOW DID WE FINISH ? */
    if( r_Anorm > tolerance ) {
        std::cout << "CRM process has failed.\n";
    } else { 
        std::cout << "CRM process has succeeded.\n";
    }

    recent_deviation = r_Anorm;
    
}
  
  
  
void ConjugateResidualMethod::iterationStart( 
    const FloatVector& x, const FloatVector& b, 
    FloatVector& r, FloatVector& d, FloatVector& Ar, FloatVector& Ad,
    Float& r_Anorm
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
    r_Anorm = r * Ar;
    assert( r_Anorm >= 0. );
      
    std::cout << "starting with"
              << " r-norm="    << r.norm()
              << " r-Anorm="   << Ar.norm() 
              << std::endl;
      
}




void ConjugateResidualMethod::iterationStep( 
    FloatVector& x, 
    FloatVector& r, FloatVector& d, FloatVector& Ar, FloatVector& Ad,
    Float& r_Anorm,
    FloatVector& AAd
) const {

    assert( r_Anorm >= 0. );
    
    /*  AAd = A * Ad */
    // A.apply( AAd, Ad, 1. );
    AAd = A * Ad;

    /*  alpha = r.A.r / d.AAd */
    Float norm_sq_Ad = Ad * Ad;
    Float alpha = r_Anorm / norm_sq_Ad;

    assert( r_Anorm >= 0. ); assert( norm_sq_Ad >= 0. );
    if( norm_sq_Ad < tolerance ) return;

    /*  x += alpha d */
    x += alpha * d;

    /*  r -= alpha Ad */
    r -= alpha * Ad; //= b - A * x;

    /*  Ar -= alpha AAd */
    Ar -= alpha * AAd;

    // std::cout << "deviation: " << ( Ar - A * r ).norm() << "\t";

    /*  beta = r.Ar / rho */
    Float tau = r_Anorm;
    r_Anorm = r * Ar; //r * ( A * r ); //r * Ar;
    Float beta = r_Anorm / tau;
    assert( r_Anorm >= 0. );
    
    /*  d = r + beta d */
    d = r + beta * d;

    /*  Ad = Ar + beta Ad */
    Ad = Ar + beta * Ad;
      
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
