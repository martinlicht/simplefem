
#include "pcrm.hpp"

#include "../operators/floatvector.hpp"

  
PreconditionedConjugateResidualMethod::PreconditionedConjugateResidualMethod( const LinearOperator& op, const LinearOperator& M )
: IterativeSolver(), A( op ), M( M )
{
    this->max_iteration_count = op.getdimin();
    PreconditionedConjugateResidualMethod::check();
}

PreconditionedConjugateResidualMethod::~PreconditionedConjugateResidualMethod()
{}
  
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
    
    Float r_MAMnorm = notanumber;
    
    FloatVector   r( dimension, 0. );
    FloatVector   p( dimension, 0. ); 
    
    FloatVector  Mr( dimension, 0. );
    FloatVector  Mp( dimension, 0. );
    FloatVector AMr( dimension, 0. );
    FloatVector AMp( dimension, 0. );
        
    // avoid repeated allocation of these temporary vectors 
    FloatVector  MAMp( dimension, 0. );    
    FloatVector AMAMp( dimension, 0. );
        
    std::cout << "Begin Preconditioned Conjugate Residual iteration" << std::endl;
    
    recent_iteration_count = 0;
    
    while(true)
    {
        
        /* Start / Restart PCRM process */
        if( recent_iteration_count % x.getdimension() == 0 ) {
        
            std::cout << "Begin Preconditioned Conjugate Residual iteration" << std::endl;
        
            iterationStart( x, b, r, p, Mr, Mp, AMr, AMp, r_MAMnorm );

            std::cout << "tolerance: " << tolerance << std::endl;

        }

        bool continue_condition = recent_iteration_count < max_iteration_count && r_MAMnorm > tolerance;
        
        /* Print information if it is time too */
        if( recent_iteration_count % print_modulo == 0 or not continue_condition ) {
            std::cout 
                << "#" << recent_iteration_count << "/" << max_iteration_count
                << " r-MAMsqnorm=" << r_MAMnorm
                << " r-Msqnorm="   << ( r * ( M * r ) ) 
                << std::endl;
        }

        /* If exit condition met, exit */
        if( not continue_condition ) 
            break;
            
        /* Perform iteration step */
        iterationStep( x, r, p, Mr, Mp, AMr, AMp, r_MAMnorm, MAMp, AMAMp );
        
        /* Increase iteration counter */
        recent_iteration_count++;
    }
        

    /* HOW DID WE FINISH ? */
    if( r_MAMnorm > tolerance ) {
      std::cout << "PCRM process has failed.\n";
    } else { 
      std::cout << "PCRM process has succeeded.\n";
    }

    recent_deviation = r_MAMnorm;
    
}
  
  
  
void PreconditionedConjugateResidualMethod::iterationStart( 
    const FloatVector& x, const FloatVector& b, 
    FloatVector& r, FloatVector& p, FloatVector& Mr, FloatVector& Mp, FloatVector& AMr, FloatVector& AMp,
    Float& r_MAMnorm
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
    r_MAMnorm = Mr * AMr;

    std::cout << "starting with"
              << " r-norm="    << r.norm()
              << " r-Mnorm="   << Mr.norm() 
              << " r-MAMnorm=" << r_MAMnorm
              << std::endl;
      
}


void PreconditionedConjugateResidualMethod::iterationStep( 
    FloatVector& x, 
    FloatVector& r, FloatVector& p, FloatVector& Mr, FloatVector& Mp, FloatVector& AMr, FloatVector& AMp,
    Float& r_MAMnorm,
    FloatVector& MAMp, FloatVector& AMAMp
) const {
    
    MAMp  = M *  AMp;
    AMAMp = A * MAMp;

    Float alpha = r_MAMnorm / ( AMp * MAMp );

    x += alpha * Mp;

     r  -= alpha *   AMp; // not necessary, see also below
     Mr -= alpha *  MAMp;
    AMr -= alpha * AMAMp;

    Float newr_MAMnorm = Mr * AMr;

    Float beta = newr_MAMnorm / r_MAMnorm;

      p =   r + beta *   p; // as such, completely useless
     Mp =  Mr + beta *  Mp;
    AMp = AMr + beta * AMp;

    r_MAMnorm = newr_MAMnorm;

}

  
  
  
  
  

