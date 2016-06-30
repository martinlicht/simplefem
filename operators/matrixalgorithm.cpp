
#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <cctype>

#include "matrixalgorithm.hpp"
#include "densematrix.hpp"
#include "floatvector.hpp"
#include "crm.hpp"



// DenseMatrix Inverse( const DenseMatrix& A ){
	
	// assert( A.issquare() );
	// DenseMatrix ret( A );
	// ret.zeromatrix();
	
	// ConjugateResidualMethod crm( A.transpose() * A );
	// crm.max_iteration_count = A.getdimout() * 10;
	// crm.error_tolerance = 1e-13;
	
	// std::cout << A << A.transpose() << A.transpose() * A << std::endl;
	
	// for( int c = 0; c < A.getdimout(); c++ ) {
		// FloatVector rhs(A.getdimout()), x(A.getdimout());
		// rhs.zero();
		// rhs[c] = 1.;
		// rhs = A.transpose() * rhs;
		// x.zero();
		// crm.solve( x, (const FloatVector&)rhs );
		
		// std::cout << A.transpose() * A * x << std::endl;
		// std::cout << A.transpose() * A * x - rhs << std::endl;
		
		// ret.setcolumn( c, A * x );
		
	// }
		
	// return ret;	
// }


void PolarDecomposition( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R )
{
	
	assert( A.getdimout() == Q.getdimout() && R.getdimout() == Q.getdimout() );
	assert( A.issquare() && Q.issquare() && R.issquare() );
	
	R.zeromatrix();
	
	for( int c = 0; c < A.getdimout(); c++ ) {
		
		FloatVector u = A.getcolumn(c);
		for( int j = 0; j < c; j++ ){
			R(j,c) = u * Q.getcolumn(j);
			u -= R(j,c) * Q.getcolumn(j);
		}
		R(c,c) = u*u;
		Q.setcolumn( c, u / R(c,c) );
		
	}
	
}










