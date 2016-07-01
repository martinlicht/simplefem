
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
		R(c,c) = sqrt( u*u );
		Q.setcolumn( c, u / R(c,c) );
	}
	
}


void PolarDecompositionRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t )
{
	if( t == 0 )
		return;
	if( t == 1 )
		PolarDecomposition( A, Q, R );
	else {
		DenseMatrix Qw(Q), Qv(Q);
		DenseMatrix Rw(R), Rv(R);
		PolarDecompositionRepeated( A, Qw, Rw, t-1 );
		PolarDecomposition( Qw, Qv, Rv );
		Q = Qv;
		R = Rv * Rw;
	}
}


Float TriangularDeterminant( const DenseMatrix& src )
{
	assert( src.issquare() );
	Float ret = 1.;
	for( int i = 0; i < src.getdimout(); i++ )
		ret *= src(i,i);
	return ret;
}

DenseMatrix UpperTriangularInverse( const DenseMatrix& src )
{
	assert( src.issquare() );
	DenseMatrix ret( src );
	ret.zeromatrix();
	const int D = src.getdimout();
	for( int c = D-1; c >= 0; c-- ) {
		ret(c,c) = 1. / src(c,c);
		for( int r = c-1; r >= 0; r-- ) {
			ret(r,c) = 0.;
			for( int k = c; k > r; k-- )
				ret(r,c) -= src(r,k) * ret(k,c);
			ret(r,c) /= src(r,r);
		}
	}
	return ret;
}


DenseMatrix CholeskyDecomposition( const DenseMatrix& src )
{
	assert( src.issquare() );
	DenseMatrix ret = src;
	ret.set( 0. );
	const int D = src.getdimout();
	for( int k = 0; k < D; k++ ){
		ret(k,k) = src(k,k);
		for( int j = 0; j < k; j++ )
			ret(k,k) -= ret(j,k) * ret(j,k);
		ret(k,k) = sqrt( ret(k,k) );
		for( int i = k+1; i < D; i++ ) {
			ret(k,i) = src(k,i);
			for( int j = 0; j < k; j++ )
				ret(k,i) -= ret(j,i) * ret(j,k);
			ret(k,i) /= ret(k,k);
		}
	}
	return ret;
}










