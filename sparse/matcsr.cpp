
#include <cassert>
#include <cmath>
#include <ostream>
#include <utility>
#include <vector>

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
    matrix.check();
    
    if( not matrix.is_sorted() ) {
        matrix.sortandcompressentries();
    }
    
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

    MatrixCSR::check();

}

MatrixCSR::MatrixCSR( int rows, int columns )
: MatrixCSR( SparseMatrix( rows, columns ) )
{
    MatrixCSR::check();
}




MatrixCSR::~MatrixCSR()
{
    MatrixCSR::check();
}





MatrixCSR::MatrixCSR( const MatrixCSR& mat )
: LinearOperator( mat.getdimout(), mat.getdimin() ),
  A( mat.A ),
  C( mat.C ),
  V( mat.V )
{
    LOG << "*************************************************\n";
    LOG << "*********** WARNING: DEEP COPY ******************\n";
    LOG << "***********  OF CSR MATRIX     ******************\n";
    LOG << "*************************************************\n";
    MatrixCSR::check();
}

MatrixCSR& MatrixCSR::operator=( const MatrixCSR& mat )
{
    LOG << "*************************************************\n";
    LOG << "********** WARNING: DEEP ASSIGN *****************\n";
    LOG << "**********    OF CSR MATRIX     *****************\n";
    LOG << "*************************************************\n";
    assert( getdimin() == mat.getdimin() );
    assert( getdimout() == mat.getdimout() );
    this->A = mat.A;
    this->C = mat.C;
    this->V = mat.V;
    check();
    return *this;
}

MatrixCSR::MatrixCSR( MatrixCSR&& mat )
: LinearOperator( mat.getdimout(), mat.getdimin() ),
  A( std::move(mat.A) ),
  C( std::move(mat.C) ),
  V( std::move(mat.V) )
{
    MatrixCSR::check();
}

MatrixCSR& MatrixCSR::operator=( MatrixCSR&& mat )
{
    assert( getdimin() == mat.getdimin() );
    assert( getdimout() == mat.getdimout() );
    this->A = std::move( mat.A );
    this->C = std::move( mat.C );
    this->V = std::move( mat.V );
    check();
    return *this;
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

std::string MatrixCSR::text() const
{
    std::string str_A, str_C, str_V; 

    for( int i = 0; i < A.size(); i++ ) str_A += ( std::to_string(A[i]) + " " );
    for( int i = 0; i < C.size(); i++ ) str_C += ( std::to_string(C[i]) + " " );
    for( int i = 0; i < V.size(); i++ ) str_V += ( std::to_string(V[i]) + " " );
    
    return std::string("CSRMatrix ") + std::to_string(getdimout()) + "x" + std::to_string(getdimin())
                        + "\nA: " + str_A
                        + "\nC: " + str_C
                        + "\nV: " + str_V;
                        // + "\n";
}

void MatrixCSR::printplain( std::ostream& os ) const
{
    print( os );
}

void MatrixCSR::apply( FloatVector& dest, const FloatVector& add, Float scaling ) const
{
    check();
    add.check();
    dest.check();
    
    assert( getdimin() == add.getdimension() );
    assert( getdimout() == dest.getdimension() );
    assert( &dest != &add );

    dest.zero();
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int i = 0; i < A.size()-1; i++ ){
        for( int j = A[i]; j < A[i+1] ; j++ ){
            dest[i] += scaling * V[j] * add[ C[j] ];
        }
    }

    // TODO: introduce inline assembler for the code above
    // -- uses outer loop assembler, inner loop assembler
    
}








void MatrixCSR::scale ( Float s )
{
    for( auto& v : this->V ) v *= s;
}

bool MatrixCSR::isfinite() const 
{
    for( const Float& value : V )
        if( not std::isfinite(value) ) 
            return false;
    return true;
}

FloatVector MatrixCSR::diagonal() const
{ 
    check();
    assert( getdimin() == getdimout() );
    auto ret = FloatVector( getdimin(), 0. );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int r = 0; r < getdimout(); r++ ) {
        
        ret[r] = 0.;
        
        for( int i = A[r]; i < A[r+1]; i++ )
            if( C[ i ] == r )
                ret[r] += V[ i ];
        
    }
    
    
        
    return ret;
}


int MatrixCSR::getnumberofzeroentries() const
{ 
    check();
    
    int ret = 0;
    
    for( const auto& value : V )
        if( value == 0. )
            ret++;
    
    return ret;

}










const int*   MatrixCSR::getA() const 
{ 
    return A.data();
}

const int*   MatrixCSR::getC() const 
{ 
    return C.data();
}

const Float* MatrixCSR::getV() const
{ 
    return V.data();
}

        
int MatrixCSR::getnumberofentries() const 
{
    return SIZECAST( V.size() );
}


// void MatrixCSR::sortentries() const
// {
//     check();
// 
//     // for each row, perform a quick sort of the columns
//     // the following uses an inefficient version of selection sort
//     for( int r = 0; r < A.size()-1; r++ )
//         for( int c1 = A[r]; c1 < A[r+1]; c1++ )
//         for( int c2 = c1+1; c2 < A[r+1]; c2++ )
//             if( C[c1] < C[c2] ) {
//                 std::swap( C[c1], C[c2] );
//                 std::swap( V[c1], V[c2] );
//             }
//     check();
// }

Float MatrixCSR::eigenvalueupperbound() const 
{
    Float ret = 0.;
    for( int r = 0; r < A.size()-1; r++ ) {
        Float candidate_ret = 0.;
        for( int c = A[r]; c < A[r+1]; c++ )
            candidate_ret += absolute( V[c] );
        ret = std::max( ret, candidate_ret );
    }
    return ret;
}







MatrixCSR MatrixCSRMultiplication( const MatrixCSR& mat1, const MatrixCSR& mat2 )
{
    // gather relevant data
    int mat1_rows = mat1.getdimout();
    int mat1_cols = mat1.getdimin();
    int mat2_rows = mat2.getdimout();
    int mat2_cols = mat2.getdimin();
    
    const int* mat1A = mat1.getA();
    const int* mat2A = mat2.getA();
    
    const int* mat1C = mat1.getC();
    const int* mat2C = mat2.getC();

    const Float* mat1V = mat1.getV();
    const Float* mat2V = mat2.getV();

    Assert( mat1_cols == mat2_rows );

    int matn_rows = mat1_rows;
    int matn_cols = mat2_cols;

    // create index list A, allocate C and V 
    
    std::vector<int> A( mat1_rows + 1, 0 );

    // LOG << mat1.text() << nl << mat2.text() << nl;

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int r = 0; r < matn_rows; r++ )
        for( int c1 = mat1A[r];           c1 < mat1A[r+1];           c1++ )
        for( int c2 = mat2A[ mat1C[c1] ]; c2 < mat2A[ mat1C[c1]+1 ]; c2++ )
            A[r+1]++;

    // for( int r = 0; r <= matn_rows; r++ ) LOG << A[r] << space; LOG << nl;
    
    for( int r = 1; r <= matn_rows; r++ ) A[r] += A[r-1];

    // for( int r = 0; r <= matn_rows; r++ ) LOG << A[r] << space; LOG << nl;
    
    // LOG << A[matn_rows] << nl;

    std::vector<int>   C( A[matn_rows] );
    std::vector<Float> V( A[matn_rows] );

    // compute entries     

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int r = 0; r < matn_rows; r++ ) {

        int c_curr = A[r];
        
        for( int c1 = mat1A[r];           c1 < mat1A[r+1];           c1++ )
        for( int c2 = mat2A[ mat1C[c1] ]; c2 < mat2A[ mat1C[c1]+1 ]; c2++ )
        {
            C[ c_curr ] = mat2C[ c2 ];
            V[ c_curr ] = mat1V[ c1 ] * mat2V[ c2 ];
            c_curr++;
        }
            
        Assert( c_curr == A[r+1], r, c_curr, A[r+1] ); 

    }
        
    // computations done, create matrix 

    return MatrixCSR( matn_rows, matn_cols, A, C, V );

}







DiagonalOperator InverseDiagonalPreconditioner( const MatrixCSR& mat )
{ 
    mat.check();
    assert( mat.getdimin() == mat.getdimout() );

    auto diag = mat.diagonal();
    
    for( int r = 0; r < mat.getdimin(); r++ )
    {
        assert( diag[r] >= 0. );
        
        if( diag[r] > 0. ) diag[r] = 1. / diag[r];
    }
    
    return DiagonalOperator( diag );
    
//     auto ret = FloatVector( mat.getdimin(), 0. );
// 
//    #if defined(_OPENMP)
//    #pragma omp parallel for
//    #endif
//    for( int r = 0; r < mat.getdimin(); r++ ) {
//         
//         ret[r] = 0.;
//         
//         for( int c = mat.A[r]; c < mat.A[r+1]; c++ )
//             if( mat.C[ c ] == r )
//                 ret[r] += mat.V[ c ];
//     
//         assert( ret[r] >= 0. );
//         
//         if( ret[r] > 0. ) ret[r] = 1. / ret[r];
//         
//     }
//      
}



































