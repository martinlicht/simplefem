
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


MatrixCSR::~MatrixCSR()
{
    MatrixCSR::check();
}



void MatrixCSR::check() const
{
    LinearOperator::check();
    
    assert( A.size() == getdimout()+1 );
    assert( A[ getdimout() ] == V.size() );

    for( int p = 1; p < getdimout(); p++ ) assert( A[p-1] <= A[p] );
    
    assert( C.size() == getdimin() );
    assert( C.size() == V.size() );
    
    for( int i = 0; i < C.size(); i++ ) assert( 0 <= C[i] && C[i] < getdimin() && std::isfinite( V[i] ) );
}

void MatrixCSR::print( std::ostream& os ) const
{
    check();
    
    os << getdimout() << ' ' << getdimin() << ' ' << V.size() << std::endl;

    for( int i = 0; i < A.size()-1; i++ ){
        os << A[i] << " ";
    }
    os << A[ A.size()-1 ] << std::endl;
    
    for( int i = 0; i < C.size()-1; i++ ){
        os << C[i] << " ";
    }
    os << C[ C.size()-1 ] << std::endl;
    
    for( int i = 0; i < V.size()-1; i++ ){
        os << V[i] << " ";
    }
    os << V[ V.size()-1 ] << std::endl;
    
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

void MatrixCSR::sortandcompressentries() const
{
    check();
    
    sortentries();
    
    // compress duplicate column entries
    
    check();  
}










// MatrixCSR operator&( const MatrixCSR& left, const MatrixCSR& right )
// {
//     left.sortandcompressentries();
//     right.sortandcompressentries();
//     
//     unreachable();
//         
//     ret.sortandcompressentries();
//     
//     return ret;
// 
// }








