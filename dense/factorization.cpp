
#include "densematrix.hpp"
#include "functions.hpp"
#include "simplesolver.hpp"
#include "factorization.hpp"

#include <new>

 // LR factorization, column pivot, 
 
DenseMatrix GaussJordan( DenseMatrix A )
{
    
    assert( A.is_square() );
    
    const int n = A.getdimout();
//     std::vector<int> pivotcol( n, -17 );
//     std::vector<int> colperm ( n, -17 );
//     for( int c = 0; c < n; c++ ) colperm[c] = c;

    
    DenseMatrix ret(n);
    ret.identity_matrix();
    
    // 1. eliminate lower triangular part, save coeffecients
    
    for( int i = 0; i < n; i++ ) {
        
        for( int k = 0; k < n; k++ ) { // each
            
            if( i == k ) continue; 
            
            Float coeff = - A( k, i ) / A( i, i );
            
            for( int j = i; j < n; j++ ) {
                A( k, j ) = A( k, j ) + coeff * A( i, j );
            }

            for( int j = 0; j <= i; j++ ) {
                ret( k, j ) = ret( k, j ) + coeff * ret( i, j );
            }
            
        }
        
        Float coeff = 1. / A(i,i);
        
        for( int j = 0; j <= i; j++ ) {
            ret(i,j) *= coeff;
        }
        
        for( int j = i; j < n; j++ ) {
            A(i,j) *= coeff;
        }
            
    }
    
// //     2. normalize the diagonals
//     
//     for( int i = 0; i < n; i++ ) {
//         
//         Float coeff = 1. / A(i,i);
//         
//         for( int k = 0; k < n; k++ ) {
//             A(i,k) *= coeff;
//             ret(i,k) *= coeff;
//         }
        
//         ret(i,i) = coeff;
//     }
    
    // finished!
    
    // LOG << A;
    // LOG << ret;
    
    return ret;
}

 

 
 
 
 
DenseMatrix GaussJordanInplace( DenseMatrix A, bool pivoting )
{
    
    assert( A.is_square() );
    
    const int n = A.getdimout();
    
    int* pivots = nullptr;
    if(pivoting) pivots = new (std::nothrow) int[n];
    
    for( int i = 0; i < n; i++ ) {
        
        if( pivoting ) {
            
            int c_max = i;
            for( int c = i+1; c < n; c++ )
                if( absolute(A(i,c)) > absolute(A(i,c_max)) ) 
                    c_max = c;
            
            pivots[i] = c_max;
            A.swapcolumn( c_max, i );
            
        }
        
        for( int k = 0; k < n; k++ ) { // each
            
            if( i == k ) continue; 
            
            assert( absolute(A(i,i)) != 0.0 );
            
            Float coeff = - A( k, i ) / A( i, i );
            
            for( int j = i+1; j < n; j++ ) {
                A( k, j ) = A( k, j ) + coeff * A( i, j );
            }

            A( k, i ) = coeff;
            
            for( int j = 0; j < i; j++ ) {
                A( k, j ) = A( k, j ) + coeff * A( i, j );
            }
            
        }
        
        Float coeff = 1. / A(i,i);
        
        for( int j = 0; j < i; j++ ) {
            A(i,j) *= coeff;
        }
        
        A(i,i) = coeff;
        
        for( int j = i+1; j < n; j++ ) {
            A(i,j) *= coeff;
        }
            
    }
    
    if( pivoting ) {
        for( int i = n-1; i >= 0; i-- )
//         for( int i = 0; i < n; i++ ) 
        {
//             LOG << "swap " << i << space << pivots[i] << nl;
            A.swaprow( i, pivots[i] );
        }
    }
    
    // finished!
    
    if( pivoting ) delete[] pivots;
    
    return A;
}








DenseMatrix CholeskyDecomposition( const DenseMatrix& A )
{
    return CholeskyDecompositionBanachchiewicz( A );
}


DenseMatrix CholeskyDecompositionBanachchiewicz( const DenseMatrix& A )
{
    A.check();
    assert( A.is_square() );

    DenseMatrix L = A;
    L.set( 0. );
    const int dim = A.getdimout();

    for( int r = 0; r < dim; r++ ){
        
        for( int c = 0; c < r; c++ ) {
        
            L(r,c) = A(r,c);
            
            for( int k = 0; k < c; k++ )
                L(r,c) -= L(r,k) * L(c,k);
            
            L(r,c) /= L(c,c);
            
            assert( absolute( L(c,c) ) > 0. );
        }
        
        L(r,r) = A(r,r);
        
        for( int k = 0; k < r; k++ )
            L(r,r) -= L(r,k) * L(r,k);
        
        assert( L(r,r) >= 0. );
        
        L(r,r) = std::sqrt( L(r,r) );
        
    }

    L.check();
    return L;
}




void QRFactorization( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R )
{
    A.check();
    Q.check();
    R.check();
    assert( A.getdimout() == Q.getdimout() );
    assert( Q.getdimin()  == R.getdimout() );
    assert( A.getdimin()  == R.getdimin()  ); 
    assert( A.getdimin()  <= A.getdimout() );
    assert( R.is_square() );
    
    R.zero_matrix();
    
    for( int c = 0; c < A.getdimin(); c++ ) {
        FloatVector u = A.getcolumn(c);
        for( int j = 0; j < c; j++ ){
                R(j,c) = u * Q.getcolumn(j);
                u -= R(j,c) * Q.getcolumn(j);
        }
        R(c,c) = std::sqrt( u*u );
        Q.setcolumn( c, u / R(c,c) );
    }
    
}

void LQFactorization( const DenseMatrix& A, DenseMatrix& L, DenseMatrix& Q )
{
    A.check();
    L.check();
    Q.check();
    
    assert( A.getdimout() == L.getdimout() );
    assert( A.getdimin()  == Q.getdimin()  ); 
    assert( Q.getdimout() == L.getdimin() );
    assert( A.getdimout() <= A.getdimin() );
    assert( L.is_square() );
    
    L.zero_matrix();
    
    for( int r = 0; r < A.getdimout(); r++ ) {
        FloatVector v = A.getrow(r);
        for( int i = 0; i < r; i++ ){
                L(r,i) = v * Q.getrow(i);
                v -= L(r,i) * Q.getrow(i);
        }
        L(r,r) = std::sqrt( v*v );
        Q.setrow( r, v / L(r,r) );
    }
    
}









// A is stored in row-major form as 1D vector of length n*n.
// That is, A[row, col] is A[row * n + col].
Float QRFactorization_via_Householder(const DenseMatrix& A_in, DenseMatrix& Q_out, DenseMatrix& R_out )
{
    assert( A_in.is_square() );

    const int n = A_in.getdimin();
    
    // Make copies so we don't destroy the input matrix.
    // R will be the in-place matrix we transform to upper triangular.
    R_out = A_in;
    
    // Q_out starts as the identity
    Q_out = IdentityMatrix(n);
    
    // Track number of nontrivial reflectors
    size_t num_reflectors = 0;

    // Householder transformations
    for( size_t k = 0; k < n; k++ )
    {
        // 1) Extract the subvector x = R[k..n-1, k]
        //    We'll store it in a local std::vector<Float> x.
        FloatVector x(n - k, 0.0);
        for( size_t i = k; i < n; i++ ) {
            x[i - k] = R_out(i,k);  // R[i, k]
        }

        // 2) Compute alpha = -sign(x[0]) * norm(x)
        Float alpha = 0.0;
        Float x0 = x[0];
        Float norm_x = x.l2norm();
        if (norm_x > 0.0) {
            alpha = -sign(x0) * norm_x;
        } else {
            // If norm_x == 0, reflection is trivial
            alpha = 0.0;
        }

        // 3) Form vector v = x; v[0] = x[0] - alpha
        //    Then normalize v so that its norm = 1 (unless it's trivial).
        x[0] = x[0] - alpha;
        Float norm_v = x.norm();

        if (norm_v > machine_epsilon) // some small tolerance
        {
            // This is a nontrivial reflector
            num_reflectors++;

            // Normalize v
            x /= norm_v;
        }
        else
        {
            // Reflection is trivial -> skip
            continue;
        }

        // 4) Apply the reflector to R on the left:
        //
        //    For each column j from k to n-1:
        //       Let w = dotProduct(v, R[k..n-1, j])
        //       Then R[k..n-1, j] -= 2 * v * w
        //
        //    Here v = x, which has length n-k, so we must offset row indices.

        for( size_t j = k; j < n; j++ )
        {
            // Compute w = dot(v, col_j_subvector)
            Float w = 0.0;
            for( size_t i = k; i < n; i++ ) {
                w += x[i - k] * R_out(i,j);  // R[i, j]
            }

            // Subtract 2 * v * w from that subvector in R
            for( size_t i = k; i < n; i++ ) {
                R_out(i,j) -= 2.0 * x[i - k] * w;
            }
        }

        // 5) If we want to build Q explicitly, 
        //    Q = Q * P_k  (where P_k = I - 2 v v^T in rows k..n-1)
        //    We need to apply the same transformation to Q_out rows [k..n-1].

        //    Recall: Q_out is n x n, v is length (n-k).
        //            For each column j in [0..n-1], 
        //            we do: Q_out[k..n-1, j] -= 2 * v * dot(v, Q_out[k..n-1, j])

        for( size_t j = 0; j < n; j++ )
        {
            Float w = 0.0;
            for( size_t i = k; i < n; i++ ) {
                w += x[i - k] * Q_out(i,j);  // Q[i, j]
            }
            for( size_t i = k; i < n; i++ ) {
                Q_out(i,j) -= 2.0 * x[i - k] * w;
            }
        }
    }

    // At this point, R_out is upper triangular, and Q_out is the orthonormal factor.

    // If you want Q explicitly as the factor such that A = Q * R,
    // the Q_out we have is indeed that matrix (since we "left-applied" Householders to R
    // while right-multiplying Q_out). 
    // In many references, you might see that Q is the product of reflectors in reverse order,
    // but here we applied them in forward order to Q_out initially set to identity.
    // Double-check consistent usage in your code base.

    // Compute determinant if desired:
    //  det(A) = det(Q) * det(R).
    //  det(Q) = (-1)^(num_reflectors) because each nontrivial Householder has det = -1.
    //  det(R) = product of diagonal entries.

    TransposeSquareInSitu(Q_out);

    // 1) diagonal product of R
    Float diag_product = 1.0;
    for( size_t i = 0; i < n; i++ ) {
        diag_product *= R_out(i,i);  // R[i, i]
    }

    // 2) sign factor from the reflectors
    Float sign_factor = ( (num_reflectors % 2) == 0 ) ? 1.0 : -1.0;

    // 3) final determinant
    return sign_factor * diag_product;
}









FloatVector SolveOverconstrained( const DenseMatrix& A, const FloatVector& b )
{
    assert( A.getdimout() == b.getdimension() );

    DenseMatrix Q( A.getdimout(), A.getdimin() );
    DenseMatrix R( A.getdimin() );
    QRFactorization(A,Q,R);

    FloatVector x( A.getdimin() );

    UpperTriangularSolve( R, x, Transpose(Q) * b );

    return x; 
}


FloatVector QRIteration( DenseMatrix A, int repetitions ) 
{
    assert( A.is_square() && A.is_symmetric() );
    
    const int dim = A.getdimin();
    
    while( repetitions --> 0 )
    {
        DenseMatrix Q(dim,dim), R(dim,dim);        
        QRFactorization( A, Q, R ); A = R * Q;
    }

    return A.getDiagonal();
}







DenseMatrix Inverse_via_LQ( const DenseMatrix& A )
{
    assert( A.is_square() );
    const int n = A.getdimin();

    // 1. Factor A = L * Q
    DenseMatrix L(n), Q(n);
    LQFactorization(A, L, Q);

    // 2. We have Ainv = Qinv * Linv = Qt * Linv 

    // 2.1. Invert the right-upper triangular matrix R
    DenseMatrix& Linv = L;
    InvertLowerTriangular(Linv);

    // 2.2. Invert the matrix Q
    DenseMatrix& Qt = Q;
    TransposeSquareInSitu(Qt);

    // 3. Return AInv = Q^t * L^-1
    return Qt * Linv;
}
DenseMatrix Inverse_via_QR( const DenseMatrix& A )
{
    assert( A.is_square() );
    const int n = A.getdimin();

    // 1. Factor A = Q * R
    DenseMatrix Q(n), R(n);
    QRFactorization(A, Q, R);

    // 2. We have Ainv = Rinv * Qinv = Rinv * Qt

    // 2.1. Invert the right-upper triangular matrix R
    DenseMatrix& Rinv = R;
    InvertUpperTriangular(Rinv);

    // 2.2. Invert the matrix Q
    DenseMatrix& Qt = Q;
    TransposeSquareInSitu(Qt);

    // 3. Return AInv = R^-1 * Qt
    return Rinv * Qt;
}


DenseMatrix Inverse_via_Cholesky(const DenseMatrix& A)
{
    // 1. Compute Cholesky factor: A = L * L^T, where L is lower-triangular.
    DenseMatrix L = CholeskyDecomposition(A);

    // 2. Invert the lower-triangular L in place.
    DenseMatrix& Linv = L;
    InvertLowerTriangular(Linv);

    // 3. Form L^-1^T by transposing the inverted L
    DenseMatrix Linvt = Transpose(Linv);

    // 4. Compute A^-1 = (L^-1)^T * L^-1
    return Linvt * Linv;
}
