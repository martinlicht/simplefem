
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
	
	/* Build up data */
    int iter = 0;
    
    Float r_Anorm;
    FloatVector& r = residual;
	FloatVector  d( dimension );  d.zero();
	FloatVector Ad( dimension ); Ad.zero();
	FloatVector Ar( dimension ); Ar.zero();
	    
    /* Initiate CRM process */
    iterationStart( x, b, residual, d, Ar, Ad, r_Anorm );
    
    std::cout << "Begin iteration" << std::endl;
	
	/* Perform CRM step */
    while( iter < max_iteration_count && r_Anorm > error_tolerance ) {
    
	  std::cout 
		<< "iteration: " << iter << "/" << max_iteration_count
		<< " : "
		<< r_Anorm << " vs " << error_tolerance
		<< std::endl;
    
      iterationStep( x, r, d, Ar, Ad, r_Anorm );
	  
      iter++;
    
	  std::cout << "iteration: " << iter << " : " << r_Anorm << std::endl;
      
    }
    
	/* FINISHED */
    if( r_Anorm > error_tolerance )
		std::cout << "CRM process has failed" << std::endl;
	recent_iteration_count = iter;
    recent_error = r_Anorm;
    
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
	Float& r_Anorm 
) const {
	
	/*  p = A * Ad */
    FloatVector  p( dimension );
	p.zero();
	internalOperator.apply( p, Ad, 1. );
      
	/*  alpha = r.A.r / d.p */
	Float alpha = r_Anorm / ( d * p );

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

