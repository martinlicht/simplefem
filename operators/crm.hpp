
#ifndef INCLUDEGUARD_CONJUGATERESIDUAL_METHOD
#define INCLUDEGUARD_CONJUGATERESIDUAL_METHOD

#include <iostream>


#include "../basic.hpp"

#include "iterativesolver.hpp"

  class ConjugateResidualMethod
  : public IterativeSolver
  {
    
    public:
	
	virtual void check() const;
	
	virtual void solve( FloatVector&, const FloatVector& ) const;
	
	ConjugateResidualMethod( const LinearOperator& op );
	virtual ~ConjugateResidualMethod();
    
	private:
	
		int dimension;
    
  };
  
  
  
  
#endif
  
  
  
  