#ifndef INCLUDEGUARD_SOLVER_ITERATIVESOLVER
#define INCLUDEGUARD_SOLVER_ITERATIVESOLVER


#include <iostream>

#include "../basic.hpp"
#include "../operators/linearoperator.hpp"


/************************
****
****  Abstract class for iterative solvers  
****  - uses iteration counter, error tolerance, and internal residual vector 
****  
************************/


  
struct IterativeSolver
{

        IterativeSolver( Float tolerance = 1.E-10, int max_iteration_count = 10, int print_modulo = 1 )
        : tolerance( tolerance ), 
          recent_deviation( 0. ), 
          max_iteration_count( max_iteration_count ),
          recent_iteration_count(0),
          print_modulo( print_modulo ) 
        {
            IterativeSolver::check();
        }

        virtual void check() const
        {
            assert( std::isfinite( tolerance ) && tolerance >= 0. );
            assert( std::isfinite( recent_deviation ) && recent_deviation >= 0. );
            
            assert( max_iteration_count >= 0 );
            assert( recent_iteration_count >= 0 );
            assert( recent_iteration_count <= max_iteration_count );

            assert( print_modulo >= 1 );
        }

        virtual void print( std::ostream& os ) const
        {
            os << "Print Iterative Solver." << std::endl;
        }

        virtual void solve( FloatVector& unknown, const FloatVector& rhs ) const = 0;

        mutable Float tolerance;
        mutable Float recent_deviation;
        
        mutable int   max_iteration_count;
        mutable int   recent_iteration_count;
        
        mutable int   print_modulo;

  };
  
 
  
  
  
  
  
#endif
  
  
  
  
