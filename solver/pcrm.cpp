
#include "pcrm.hpp"

#include "../operators/floatvector.hpp"

  
void PreconditionedConjugateResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( internalOperator.getdimin()  == x.getdimension() );
    assert( internalOperator.getdimout() == b.getdimension() );
    
    /* Build up data */
    int iter = 0;
    
    Float r_MAMnorm = notanumber;
    
    FloatVector&  r = residual;
    FloatVector   p( dimension, 0. ); 
    
    FloatVector  Mr( dimension, 0. );
    FloatVector  Mp( dimension, 0. );
    FloatVector AMr( dimension, 0. );
    FloatVector AMp( dimension, 0. );
        
    // avoid repeated allocation of these temporary vectors 
    FloatVector  MAMp( dimension, 0. );    
    FloatVector AMAMp( dimension, 0. );
        
    iterationStart( x, b, r, p, Mr, Mp, AMr, AMp, r_MAMnorm ); /* Initiate PCRM process */
    
    std::cout << "Begin Preconditioned Conjugate Residual iteration" << std::endl;
    std::cout << "start: " << r_MAMnorm << " tolerance: " << tolerance << std::endl;

    while( iter < max_iteration_count && r_MAMnorm > tolerance ) /* Perform PCRM step */
    {
        if( iter % print_modulo == 0 ) 
          std::cout 
            << "#" << iter << "/" << max_iteration_count
            << " : "
            << r_MAMnorm
            << " : "
            << ( r * ( invprecon * r ) ) 
            << std::endl;
            
        iterationStep( x, b, r, p, Mr, Mp, AMr, AMp, r_MAMnorm, MAMp, AMAMp );
        
        iter++;
    }
    
    /* HOW DID WE FINISH ? */
    if( r_MAMnorm > tolerance ) {
      std::cout << "PCRM process has failed. ";
    } else { 
      std::cout << "PCRM process has succeeded. ";
    }

    std::cout
       << "iterations "
       << iter << "/" << max_iteration_count 
       << " : " << r_MAMnorm << " vs " << tolerance
       << std::endl;
    
        
    iterationStart( x, b, r, p, Mr, Mp, AMr, AMp, r_MAMnorm );
    
    std::cout 
        << "recomputed "
        << r_MAMnorm << " vs " << tolerance
        << std::endl;
        
    recent_iteration_count = iter;
    recent_deviation = r_MAMnorm;
    
}
  
  
  
void PreconditionedConjugateResidualMethod::iterationStart( 
    const FloatVector& x, const FloatVector& b, 
    FloatVector& r, FloatVector& p, FloatVector& Mr, FloatVector& Mp, FloatVector& AMr, FloatVector& AMp,
    Float& r_MAMnorm
) const {
    
    const LinearOperator& A = internalOperator;
    const LinearOperator& M = invprecon;

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

    std::cout << "starting with " << r.norm() << space << Mr.norm() << space << r_MAMnorm << std::endl;
      
}


void PreconditionedConjugateResidualMethod::iterationStep( 
    FloatVector& x, const FloatVector& b,
    FloatVector& r, FloatVector& p, FloatVector& Mr, FloatVector& Mp, FloatVector& AMr, FloatVector& AMp,
    Float& r_MAMnorm,
    FloatVector& MAMp, FloatVector& AMAMp
) const {
    
    const LinearOperator& A = internalOperator;
    const LinearOperator& M = invprecon;

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

  
  
  
  
  

PreconditionedConjugateResidualMethod::PreconditionedConjugateResidualMethod( const LinearOperator& op, const LinearOperator& invprecon )
: IterativeSolver( op ), dimension( op.getdimout() ),
  invprecon( invprecon )
{
    PreconditionedConjugateResidualMethod::check();
}

PreconditionedConjugateResidualMethod::~PreconditionedConjugateResidualMethod()
{
	
}
    
	
  
void PreconditionedConjugateResidualMethod::check() const
{
    const LinearOperator& op = internalOperator;
    
    IterativeSolver::check();
    assert( op.getdimin() == op.getdimout() );
    assert( invprecon.getdimin() == invprecon.getdimout() );
    assert( op.getdimin() == invprecon.getdimout() );
}

void PreconditionedConjugateResidualMethod::print( std::ostream& os ) const
{
    os << "Print Preconditioned Conjugate Residual Methop." << std::endl;
}

