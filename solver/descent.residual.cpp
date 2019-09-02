
#include "descent.residual.hpp"

#include "../operators/floatvector.hpp"


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

 


  
void ResidualDescentMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == dimension );
    assert( b.getdimension() == dimension );
    
    /* Build up data */
    int iter = 0;
    
    Float r_norm;
    FloatVector& r = residual;
    
    /* Initiate Residual Descent process */
    iterationStart( x, b, r, r_norm );
    
    std::clog << "Begin iteration" << std::endl;

    /* Perform Residual Descent step */
    while( iter < max_iteration_count && r_norm > error_tolerance )
    {
        std::clog 
          << "iteration: " << iter << "/" << max_iteration_count
          << " : "
          << r_norm << " vs " << error_tolerance
          << std::endl;
            
        iterationStep( x, r, r_norm );
        
        FloatVector test = internalOperator * x - b;
        Float r_norm_test = test * test;
        std::clog << "[ratio: " << r_norm / r_norm_test << "]" << std::endl;
        
        iter++;
    }
    
    /* FINISHED */
    if( r_norm > error_tolerance ) {
      
      std::clog << "CRM process has failed" << std::endl;
      
    } else {
        
      iterationStart( x, b, residual, r_norm );
      std::clog 
        << "iteration: " << iter << "/" << max_iteration_count
        << " : "
        << r_norm << " vs " << error_tolerance
        << std::endl;
        
    }
    
    recent_iteration_count = iter;
    recent_error = r_norm;
    
}
  
  
  
void ResidualDescentMethod::iterationStart( 
    const FloatVector& x, const FloatVector& b, 
    FloatVector& r, 
    Float& r_norm
) const {
    
    /* r = b - A x */
    r = b - internalOperator * x;

    /* rho is r.A.r */
    r_norm = r * r;
      
}


void ResidualDescentMethod::iterationStep( 
    FloatVector& x,
    FloatVector& r, 
    Float& r_norm
) const {
    
    /*  Ar = A * r */
    FloatVector Ar = internalOperator * r;
  
    /*  alpha = r.A.r / d.p */
    Float alpha = r * Ar / ( Ar * Ar );

    /*  x += alpha d */
    x += alpha * r;

    /*  r -= alpha Ad */
    r -= alpha * Ar;
    
    r_norm = r*r;
      
}

  
  
  
  
  
  
ResidualDescentMethod::ResidualDescentMethod( const LinearOperator& op )
: IterativeSolver( op ), dimension( op.getdimout() )
{
    assert( op.getdimin() == op.getdimout() );
    ResidualDescentMethod::check();
}

ResidualDescentMethod::~ResidualDescentMethod()
{
        
}
    
        
  
void ResidualDescentMethod::check() const
{
    IterativeSolver::check();
    assert( getdimin() == getdimout() );
    assert( getdimin() == dimension );
}

void ResidualDescentMethod::print( std::ostream& os ) const
{
    os << "Print Conjugate Residual Method." << std::endl;
}


