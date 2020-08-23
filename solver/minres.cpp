
#include "minres.hpp"

#include "../operators/floatvector.hpp"



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
    
    while(true)
    {
        
        /* Start / Restart MinimumResidualMethod process */
        if( recent_iteration_count % x.getdimension() == 0 ) {
        
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
        
            if( verbosity >= VerbosityLevel::verbose ) 
            LOG << "starting with"
                      << " r-sqnorm="    << rr 
                      ;//<< std::endl;
            if( verbosity >= VerbosityLevel::verbose ) 
            LOG << "tolerance: " << tolerance;// << std::endl;

        }

        bool continue_condition = recent_iteration_count < max_iteration_count && rr > tolerance;
        
        /* Print information if it is time too */
        if( verbosity >= VerbosityLevel::verbose ) 
        if( recent_iteration_count % print_modulo == 0 or not continue_condition ) {
            LOG 
                << "#" << recent_iteration_count << "/" << max_iteration_count
                << " r-sqnorm="  << rr 
                ;//<< std::endl;
        }

        /* If exit condition met, exit */
        if( not continue_condition ) 
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
    recent_deviation = rr;
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
