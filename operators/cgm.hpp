
#ifndef CLASS_CONJUGATEGRADIENT_SOLVER
#define CLASS_CONJUGATEGRADIENT_SOLVER

#include <iostream>

using namespace std;

#include "../basic/basic.hpp"

#include "vectorspace.hpp"
#include "linearoperator.hpp"

#include "iterativesolver.hpp"

  class conjugategradientsolver
  : public iterativesolver
  {
    
  public:
    
    conjugategradientsolver( const vectorspace* V, const linearoperator* A )
    : iterativesolver( V, A )
    {
      
    }
    
    virtual ~conjugategradientsolver()
    {
      
    }
    
    virtual void run() const override;
    
  };
  
  
  
  
#endif
  
  
  
  