
#include "iterativesolver.hpp"

#include "scalingoperator.hpp"
	
IterativeSolver::IterativeSolver( const LinearOperator& op )
: LinearOperator( op.getdimout(), op.getdimin() ), 
  internalOperator( op ), 
  residual(op.getdimout()), 
  error_tolerance( 1.E-10 ), 
  recent_error( 0. ), 
  max_iteration_count(100),
  recent_iteration_count(0)
{
  check();
}

IterativeSolver::~IterativeSolver()
{}

void IterativeSolver::check() const
{
  LinearOperator::check();
  
  assert( error_tolerance >= 0. );
  assert( max_iteration_count >= 0 );
  assert( recent_iteration_count >= 0 );
  assert( recent_iteration_count <= max_iteration_count );
  
  assert( internalOperator.getdimout() == getdimout() );
  assert( internalOperator.getdimin()  == getdimin() );
  
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


void IterativeSolver::applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const
{
	FloatVector temp( getdimout() );
	temp.zero();
	solve( temp, add );
	ScalingOperator unitmatrix( getdimout(), 1. );
	unitmatrix.applyadd( dest, temp, s, t );
}
	
	

    