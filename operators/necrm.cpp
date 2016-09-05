
#include "crm.hpp"
#include "necrm.hpp"

#include "floatvector.hpp"


NormalEquationsConjugateResidualMethod::NormalEquationsConjugateResidualMethod( const DenseMatrix& op )
: matrix_original( op ), 
  matrix_transposed( matrix_original.transpose() ), 
  matrix_system( matrix_transposed * matrix_original ),
  crm( matrix_system )
{
    
}

NormalEquationsConjugateResidualMethod::~NormalEquationsConjugateResidualMethod()
{
    /* TODO */
}

  
void NormalEquationsConjugateResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    
    check();
    assert( x.getdimension() == matrix_original.getdimin() );
    assert( b.getdimension() == matrix_original.getdimout() );
    
    crm.solve( x, matrix_transposed * b );
    
}
  
  

  
	
  
void NormalEquationsConjugateResidualMethod::check() const
{
	
}

void NormalEquationsConjugateResidualMethod::print( std::ostream& os ) const
{
	os << "Print Normal Equations Conjugate Residual Method." << std::endl;
}


