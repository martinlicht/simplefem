
#include "cgm.hpp"

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
    
    while(true)
    {
        
        /* Start / Restart CRM process */
        if( recent_iteration_count % x.getdimension() == 0 ) {
        
            LOG << "Begin Conjugate Gradient iteration";// << std::endl;
        
            r = b - A * x;
            d = A * r;
            
            LOG << "starting with"
                      << " r-sqnorm="    << r * r 
                      ;//<< std::endl;
            LOG << "tolerance: " << tolerance;// << std::endl;

        }

        bool continue_condition = recent_iteration_count < max_iteration_count && r * r > tolerance;
        
        /* Print information if it is time too */
        if( recent_iteration_count % print_modulo == 0 or not continue_condition ) {
            LOG 
                << "#" << recent_iteration_count << "/" << max_iteration_count
                << " r-sqnorm="  << r * r 
                ;//<< std::endl;
        }

        /* If exit condition met, exit */
        if( not continue_condition ) 
            break;
            
        /* Perform iteration step */
        {
            
            Ad = A * d;
            
            Float rr_old = r * r; assert( rr_old >= 0 ); if( rr_old < tolerance ) break;
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
    recent_deviation = r * r;
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
