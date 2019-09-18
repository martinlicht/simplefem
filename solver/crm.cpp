
#include "crm.hpp"

#include "../operators/floatvector.hpp"

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
  

  
void ConjugateResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == dimension );
    assert( b.getdimension() == dimension );
    
    /* Build up data */
    int iter = 0;
    
    Float r_Anorm = notanumber;

    FloatVector&  r = residual;
    FloatVector   d( dimension, 0. );
    
    FloatVector  Ad( dimension, 0. );
    FloatVector  Ar( dimension, 0. );
    
    // avoid repeated allocation of these temporary vectors 
    FloatVector AAd( dimension, 0. );
    
    iterationStart( x, b, residual, d, Ar, Ad, r_Anorm ); /* Initiate CRM process */
    
    std::cout << "Begin Conjugate Residual iteration" << std::endl;
    std::cout << "start: " << r_Anorm << "tolerance: " << tolerance << std::endl;

    while( iter < max_iteration_count && r_Anorm > tolerance ) /* Perform CRM step */
    {
        if( iter % print_modulo == 0 ) 
          std::cout 
            << "#" << iter << "/" << max_iteration_count
            << " : "
            << r_Anorm 
            << std::endl;
            
        iterationStep( x, r, d, Ar, Ad, r_Anorm, AAd );
        
        iter++;
    }
    
    /* HOW DID WE FINISH ? */
    if( r_Anorm > tolerance ) {
      std::cout << "CRM process has failed. ";
    } else { 
      std::cout << "CRM process has succeeded. ";
    }

    std::cout
       << "iterations "
       << iter << "/" << max_iteration_count 
       << " : " << r_Anorm << " vs " << tolerance
       << std::endl;
    
    
    iterationStart( x, b, residual, d, Ar, Ad, r_Anorm );
    std::cout 
        << "iteration: " << iter << "/" << max_iteration_count
        << " : "
        << r_Anorm << " vs " << tolerance
        << std::endl;
    
    recent_iteration_count = iter;
    recent_deviation = r_Anorm;
    
}
  
  
  
void ConjugateResidualMethod::iterationStart( 
    const FloatVector& x, const FloatVector& b, 
    FloatVector& r, FloatVector& d, FloatVector& Ar, FloatVector& Ad,
    Float& r_Anorm
) const {
    
    /* r = b - A x */
    r = b - internalOperator * x;

    /* d = r */
    d.copydatafrom( r );

    /* Ar = A r */
    internalOperator.apply( Ar, (const FloatVector&) r );

    /* Ad = A d */
    internalOperator.apply( Ad, (const FloatVector&) d );

    /* rho is r.A.r */
    r_Anorm = r * Ar;
      
}


void ConjugateResidualMethod::iterationStep( 
    FloatVector& x,
    FloatVector& r, FloatVector& d, FloatVector& Ar, FloatVector& Ad,
    Float& r_Anorm,
    FloatVector& p
) const {
    
    /*  p = A * Ad */
    internalOperator.apply( p, Ad, 1. );
  
    /*  alpha = r.A.r / d.p */
    Float alpha = r_Anorm / ( Ad * Ad );

    /*  x += alpha d */
    x += alpha * d;

    /*  r -= alpha Ad */
    r -= alpha * Ad;

    /*  Ar -= alpha p */
    Ar -= alpha * p;

    /*  beta = r.Ar / rho */
    Float tau = r_Anorm;
    r_Anorm = r * Ar;
    Float beta = r_Anorm / tau;

    /*  d = r + beta d */
    d = r + beta * d;

    /*  Ad = Ar + beta Ad */
    Ad = Ar + beta * Ad;
      
}

  
  
  
  
  
  
ConjugateResidualMethod::ConjugateResidualMethod( const LinearOperator& op )
: IterativeSolver( op ), dimension( op.getdimout() )
{
    ConjugateResidualMethod::check();
}

ConjugateResidualMethod::~ConjugateResidualMethod()
{
	
}
    
	
  
void ConjugateResidualMethod::check() const
{
    const LinearOperator& op = internalOperator;
    
    IterativeSolver::check();
    assert( op.getdimin() == op.getdimout() );
    assert( getdimin() == op.getdimin() );
    assert( getdimin() == getdimout() );
    assert( getdimin() == dimension );
}

void ConjugateResidualMethod::print( std::ostream& os ) const
{
    os << "Print Conjugate Residual Method." << std::endl;
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
