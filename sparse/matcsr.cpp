
#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <vector>
#include <cctype>

#include "matcsr.hpp"



MatrixCSR::MatrixCSR( 
    int rows,
    int columns,
    const std::vector<int>& A, 
    const std::vector<int>& C, 
    const std::vector<Float>& V
): LinearOperator( rows, columns ),
   A(A), C(C), V(V) 
{
    MatrixCSR::check();
}


MatrixCSR::MatrixCSR( 
    const SparseMatrix& matrix
): LinearOperator( matrix.getdimout(), matrix.getdimin() ),
   A(0), C(0), V(0) 
{
    assert( matrix.is_sorted() );
    
    int rows       = matrix.getdimout();
    int columns    = matrix.getdimin();
    int numentries = matrix.getnumberofentries();

    // std::vector<int>   A( rows+1,     0  );
    // std::vector<int>   C( numentries, 0  );
    // std::vector<Float> V( numentries, 0. );
    
    A.resize( rows+1     );
    C.resize( numentries );
    V.resize( numentries );
    
    for( int i = 0; i < matrix.getnumberofentries(); i++ ){
        A[ matrix.getentry(i).row+1 ] += 1;
        C[i] = matrix.getentry(i).column;
        V[i] = matrix.getentry(i).value;
    }
 
    for( int i = 1; i < A.size(); i++ ){
        A[i] += A[i-1];
    }

    check();

}


MatrixCSR::~MatrixCSR()
{
    MatrixCSR::check();
}



void MatrixCSR::check() const
{
    LinearOperator::check();
    
    // check array dimensions and some fix values
    assert( A.size() == getdimout()+1 );
    assert( C.size() == V.size() );
    
    // check that A is ascending 
    for( int p = 1; p <= getdimout(); p++ ) assert( A[p-1] <= A[p] );
    
    // chekc the final values of A and the validity of the values in C and V
    assert( A[ getdimout() ] == V.size() );
    assert( A[ getdimout() ] == C.size() );
    for( int i = 0; i < C.size(); i++ ) assert( 0 <= C[i] && C[i] < getdimin() && std::isfinite( V[i] ) );
}

void MatrixCSR::print( std::ostream& os ) const
{
    check();
    
    os << getdimout() << ' ' << getdimin() << ' ' << V.size() << std::endl;

    for( int i = 0; i < A.size(); i++ ) os << A[i] << " ";
    os << std::endl;
    
    for( int i = 0; i < C.size(); i++ ) os << C[i] << " ";
    os << std::endl;
    
    for( int i = 0; i < V.size(); i++ ) os << V[i] << " ";
    os << std::endl;
    
}

void MatrixCSR::printplain( std::ostream& os ) const
{
    print( os );
}

FloatVector MatrixCSR::apply( const FloatVector& add, Float scaling ) const
{
    check();
    add.check();
    assert( getdimin() == add.getdimension() );

    FloatVector ret( getdimout(), 0 );  

    for( int i = 0; i < A.size()-1; i++ ){
        for( int j = A[i]; j < A[i+1] ; j++ ){
            ret[i] += scaling * V[j] * add[ C[j] ];
        }
    }
    return ret;
}







int MatrixCSR::getnumberofentries() const 
{
    return V.size();
}


void MatrixCSR::sortentries() const
{
    check();

    // for each row, perform a quick sort of the columns
    // the following uses an inefficient version of selection sort
    for( int r = 0; r < A.size()-1; r++ )
        for( int c1 = A[r]; c1 < A[r+1]; c1++ )
        for( int c2 = c1+1; c2 < A[r+1]; c2++ )
            if( C[c1] < C[c2] ) {
                std::swap( C[c1], C[c2] );
                std::swap( V[c1], V[c2] );
            }
    check();
}


















