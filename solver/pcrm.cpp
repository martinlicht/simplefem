
#include "pcrm.hpp"

#include "../operators/floatvector.hpp"

  
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
        
            iterationStart( x, b, r, p, Mr, Mp, AMr, AMp, rMAMr );

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
        iterationStep( x, r, p, Mr, Mp, AMr, AMp, rMAMr, MAMp, AMAMp );
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
        

    /* HOW DID WE FINISH ? */
    if( rMAMr > tolerance ) {
      LOG << "PCRM process has failed.\n";
    } else { 
      LOG << "PCRM process has succeeded.\n";
    }

    recent_deviation = rMAMr;
    
}
  
  
  
void PreconditionedConjugateResidualMethod::iterationStart( 
    const FloatVector& x, const FloatVector& b, 
    FloatVector& r, FloatVector& p, FloatVector& Mr, FloatVector& Mp, FloatVector& AMr, FloatVector& AMp,
    Float& rMAMr
) const {
    
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


void PreconditionedConjugateResidualMethod::iterationStep( 
    FloatVector& x, 
    FloatVector& r, FloatVector& p, FloatVector& Mr, FloatVector& Mp, FloatVector& AMr, FloatVector& AMp,
    Float& rMAMr,
    FloatVector& MAMp, FloatVector& AMAMp
) const {
    
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

  
  
  
  

