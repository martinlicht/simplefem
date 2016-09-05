
#include "crm.hpp"
#include "necrm.hpp"

#include "floatvector.hpp"


NormalEquationsConjugateResidualMethod::NormalEquationsConjugateResidualMethod( const DenseMatrix& op, const DenseMatrix& opt, const DenseMatrix& opsys )
: matrix_original( op ), 
  matrix_transposed( opt ), 
  matrix_system( opsys ),
  crm( opsys )
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


