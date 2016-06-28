
#include "crm.hpp"
  

/* 
 * 
 * Input: A, b, x
 * initialization
 * r := b - A x
 * d := r
 * Ar = A r
 * Ad = A d
 * rho = r . Ar
 * loop until exit condition
 *   if restart then reinitialize
 *   p = A * Ad
 *   alpha = rho / d . p
 *   x += alpha d
 *   r -= alpha Ad
 *   Ar -= alpha p
 *   tau = rho
 *   rho = r. Ar
 *   beta = rho / tau
 *   d = r + beta d
 *   Ad = Ar + beta Ad
 * 
 */  
  

  
  void ConjugateResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
  {
    
	check();
    assert( x.getdimension() == dimension );
	assert( b.getdimension() == dimension );
	
	int restart_period = 0;
	
	/* Build up data */
    
	int iter = 0;
    
    Float rho;
    
	FloatVector  d( dimension );  d.zero();
	FloatVector Ad( dimension ); Ad.zero();
	FloatVector Ar( dimension ); Ar.zero();
	FloatVector  p( dimension );  p.zero();
	    
    /* Algorithm */
    
    /* Initialize */
    INITIALIZE:
    {
    
      std::cout << "Initialize" << std::endl;
      
      /* r = b - A x */
      residual = b - internalOperator * x;
	  
      /* d = r */
      d.copydatafrom( residual );

      /* Ar = A r */
      internalOperator.apply( Ar, (const FloatVector&) residual );
      
      /* Ad = A d */
      internalOperator.apply( Ad, (const FloatVector&) d );
      
      /* rho is r.A.r */
      rho = residual * Ar;
    
    }
    
    /* Main iteration */
    
	std::cout << "iteration: " << iter << "/" << max_iteration_count << " : " << rho << " vs " << error_tolerance << std::endl;
      
      
	
    /* while keep running */
    while( iter < max_iteration_count && rho > error_tolerance ) {
    
      std::cout << "iteration: " << iter << " : " << rho << std::endl;
      
      /* if restart condition holds, then jump back */
      if( restart_period != 0 && iter % restart_period == 0 )
		goto INITIALIZE; // MUAHAHAHAHAHA!!!!!!!!!!1111111
      
      /*  p = A * Ad */
      internalOperator.apply( p, Ad, 1. );
      
      /*  alpha = r.A.r / d.p */
      Float alpha = rho / ( d * p );
    
      /*  x += alpha d */
      x += alpha * d;
	  
      /*  r -= alpha Ad */
      residual -= alpha * Ad;
	  
      /*  Ar -= alpha p */
      Ar -= alpha * p;
	  
	  /*  beta = r.Ar / rho */
      Float tau = rho;
      rho = residual * Ar;
      Float beta = rho / tau;
      
      /*  d = r + beta d */
      d = residual + beta * d;
	  
      /*  Ad = Ar + beta Ad */
      Ad = Ar + beta * d;
	  
      iter++;
    
	  std::cout << "iteration: " << iter << " : " << rho << std::endl;
      
      
    }
    /* FINISHED */
    
    recent_iteration_count = iter;
    recent_error = rho;
    
  }
  
  
ConjugateResidualMethod::ConjugateResidualMethod( const LinearOperator& op )
: IterativeSolver( op ), dimension( op.getdimout() )
{
	assert( op.getdimin() == op.getdimout() );
	check();
}

ConjugateResidualMethod::~ConjugateResidualMethod()
{
	
}
    
	
  
void ConjugateResidualMethod::check() const
{
	IterativeSolver::check();
	assert( getdimin() == getdimout() );
}

void ConjugateResidualMethod::print( std::ostream& os ) const
{
	os << "Print Conjugate Residual Method." << std::endl;
}

