
#include "resdes.hpp"

#include "../operators/floatvector.hpp"


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
    if( r_sqnorm > tolerance ) {
        LOG << "Residual descent process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ")\n";
    } else { 
        LOG << "Residual descent process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ")\n";

    }

    recent_deviation = r_sqnorm;
    
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

 


  
  
  
