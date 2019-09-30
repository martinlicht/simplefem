
#include <cmath>

#include "../operators/simpleoperators.hpp"
#include "iterativesolver.hpp"


IterativeSolver::IterativeSolver( const LinearOperator& op )
: internalOperator( op ), 
  residual(op.getdimout()), 
  tolerance( 1.E-10 ), 
  recent_deviation( 0. ), 
  max_iteration_count(op.getdimout()),
  recent_iteration_count(0),
  print_modulo( 1 ) 
 {
  IterativeSolver::check();
}

IterativeSolver::~IterativeSolver()
{}

void IterativeSolver::check() const
{
    assert( std::isfinite( tolerance ) && tolerance >= 0. );
    assert( std::isfinite( recent_deviation ) && recent_deviation >= 0. );
    assert( max_iteration_count >= 0 );
    assert( recent_iteration_count >= 0 );
    assert( recent_iteration_count <= max_iteration_count );
    assert( print_modulo >= 1 );
    
    internalOperator.check();
    residual.check();
}

void IterativeSolver::print( std::ostream& os ) const
{
    os << "Print Iterative Solver." << std::endl;
}

const LinearOperator& IterativeSolver::getInternalOperator() const
{
    return internalOperator;
}

const FloatVector& IterativeSolver::getResidualVector() const
{
    return residual;
}

