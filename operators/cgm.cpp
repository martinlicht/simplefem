
#include "cgm.hpp"
  
  
// \begin{algorithm}
// \caption{Conjugate Gradient}
// \begin{algorithmic}[1]
//  \REQUIRE $A \in \Reals^{N\times N}, b \in \Reals^N, x \in \Reals^N$
//  \STATE $r \leftarrow b - A \cdot x$
//  \STATE $d \leftarrow r$
//  \WHILE{ $\| r \| > \epsilon$ }
//  \STATE $t \leftarrow \| r \|$
//  \STATE $p \leftarrow A \cdot d$
//  \STATE $\alpha \leftarrow t / ( d \cdot p )$
//  \STATE $ x \leftarrow x + \alpha d$
//  \STATE $ r \leftarrow r - \alpha p$
//  \STATE $\beta = \| r \| / t$
//  \STATE $d \leftarrow r + \beta d$
//  \ENDWHILE
//  \end{algorithmic}
// \end{algorithm}
  
  
  void conjugategradientsolver::run() const
  {
    check();
    
    unsigned int iter = 0;
    const vectorspace*    space = getspace();
    const linearoperator*   ops = getlinearoperator();
    
    /* Build up data */
    Float rho;
    vectorelement* d = space->newvectorelement(); d->zero();
    vectorelement* p = space->newvectorelement(); p->zero();
    
    /* Algorithm */
    
    /* Initialize */
    INITIALIZE:
    {
    
      cout << "Initialize" << endl;
      
      /* r = b - A x */
      r->copy( b, 1.); ops->applyadd( r, x, 1., -1. );
      
      /* d = r */
      d->copy( r, 1. );
      
      /* rho is r.r */
      rho = space->scalarproduct( r, r );
    
    }
    
    /* Main iteration */
    
    /* while keep running */
    while( iter < maximal_iterations && rho > error_tolerance ) {
    
      cout << "iteration: " << iter << " : " << rho << endl;
      
      /* if restart condition holds, then jump back */
      if( restart_period != 0 && iter % restart_period == 0 )
	goto INITIALIZE; // MUAHAHAHAHAHA!!!!!!!!!!11
      
      /*  p = A d */
      ops->applyadd( p, d, 0., 1. );
      
      /*  alpha = norm of r / d . p */
      Float alpha = rho / space->scalarproduct( p, d );
    
      cout << "alpha = " << alpha << endl;
      
      /*  x += alpha d */
      x->add( d, 1.,  alpha );
      
      /*  r -= alpha p */
      r->add( p, 1., -alpha );
    
      /*  beta = r.r / rho */
      Float tau = rho;
      rho = space->scalarproduct( r, r );
      Float beta = rho / tau;
      
      /*  d = r + beta d */
      d->add( r, beta, 1. );
    
      iter++;
      
    }
    /* FINISHED */
    
    cout << "iteration: " << iter << " : " << rho << endl;
    
    /* Tear down data */
    delete d;
    delete p;
    
    last_iterations = iter;
    last_error = rho;
    
  }
  
  
  
  
  