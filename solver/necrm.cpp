
#include "crm.hpp"
#include "necrm.hpp"

#include "../operators/floatvector.hpp"


NormalEquationsConjugateResidualMethod::NormalEquationsConjugateResidualMethod( const DenseMatrix& op, const DenseMatrix& opt, const DenseMatrix& opsys )
: matrix_original( op ), 
  matrix_transposed( opt ), 
  matrix_system( opsys ),
  crm( opsys )
{
    
}

NormalEquationsConjugateResidualMethod::~NormalEquationsConjugateResidualMethod()
{
    /* Nothing to do here */
}



void NormalEquationsConjugateResidualMethod::check() const
{
  matrix_original.check();
  matrix_transposed.check();
  matrix_system.check();
  crm.check();
}

void NormalEquationsConjugateResidualMethod::print( std::ostream& os ) const
{
  os << "Print Normal Equations Conjugate Residual Method." << std::endl;
}

  
void NormalEquationsConjugateResidualMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == matrix_original.getdimin() );
    assert( b.getdimension() == matrix_original.getdimout() );
    
    crm.solve( x, matrix_transposed * b );
    
}
  
  

  


