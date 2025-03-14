
#include <cmath>

#include <functional>

#include "../../base/include.hpp"
#include "../../dense/factorization.hpp"
#include "../../dense/functions.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test: Matrix Algorithms" << nl;
    
    // LOG << std::setprecision(5);
    // LOG << std::fixed << std::ios::floatfield; //setf( std::ios::fixed, std::ios::floatfield ); 
    // LOG << std::showpos;
    
    {
        LOG << "1. Scalar functions of matrices" << nl;
       
        DenseMatrix S( 4, 4 );
        S(0,0) =  3; S(0,1) =  0; S(0,2) = -6; S(0,3) =  0; 
        S(1,0) =  1; S(1,1) = -1; S(1,2) =  0; S(1,3) =  2; 
        S(2,0) = -1; S(2,1) =  1; S(2,2) =  1; S(2,3) = -1; 
        S(3,0) =  2; S(3,1) = -4; S(3,2) =  4; S(3,3) =  0;
        
        
        LOG << S << nl;
        LOG << "Matrix trace:   " << S.trace() << nl;
        LOG << "Norm L1:        " << S.sumnorm() << nl;
        LOG << "Norm Frobenius: " << S.frobeniusnorm() << nl;
        LOG << "Norm Max:       " << S.maxnorm() << nl;

        assert( S.trace()         == 3-1+1+0 );
        assert( S.sumnorm()       == 7+6+11+3 );
        assert( S.frobeniusnorm() == std::sqrt((Float)91) );
        assert( S.maxnorm()       == 6 );
        
        Float p = 1.01;
        LOG << nl << "Norm Lp with p=" << p << ": " << S.lpnorm( p ) << nl << nl;
        
        Float p1 = 200.00000001;
        Float p2 = 1.0000000001;
        LOG << nl << "Row " << p1 << space
                     << "Col " << p2 << space
                     << S.norm_row_col( p1, p2 ) << nl;
        LOG << nl << "Col " << p1 << space
                     << "Row " << p2 << space 
                     << S.norm_col_row( p1, p2 ) << nl;
        
        Assert( is_numerically_close( 10, S.norm_row_col( p1, p2 ), 1e-7 ), S.norm_row_col( p1, p2 ) );
        Assert( is_numerically_close( 11, S.norm_col_row( p1, p2 ), 1e-7 ), S.norm_col_row( p1, p2 ) );
        
        
        LOG << nl << "Row " << 1. << space
                     << "Col " << 1. << space
                     << S.norm_row_col( 1., 1. ) << nl;
        LOG << nl << "Col " << 1. << space
                     << "Row " << 1. << space 
                     << S.norm_col_row( 1., 1. ) << nl;
        LOG << nl;
        
        assert( is_numerically_close( 27, S.norm_row_col( 1., 1. ) ) );
        assert( is_numerically_close( 27, S.norm_col_row( 1., 1. ) ) );
        
        LOG << nl << "Row " << 2. << space
                     << "Col " << 2. << space
                     << S.norm_row_col( 2., 2. ) << nl;
        LOG << nl << "Col " << 2. << space
                     << "Row " << 2. << space 
                     << S.norm_col_row( 2., 2. ) << nl;
        LOG << nl;

        Assert( is_numerically_close( S.frobeniusnorm(), S.norm_row_col( 2., 2. ) ) );
        Assert( is_numerically_close( S.frobeniusnorm(), S.norm_col_row( 2., 2. ) ) );
        
        LOG << nl << "Row " << 20. << space
                     << "Col " << 20. << space
                     << S.norm_row_col( 20., 20. ) << nl;
        LOG << nl << "Col " << 20. << space
                     << "Row " << 20. << space 
                     << S.norm_col_row( 20., 20. ) << nl;
        LOG << nl;

        Assert( is_numerically_close( 6, S.norm_row_col( 200., 200. ) ), S.norm_row_col( 200., 200. ) );
        Assert( is_numerically_close( 6, S.norm_col_row( 200., 200. ) ), S.norm_col_row( 200., 200. ) );
        
        
        
        
        
        LOG << "Norm Operator L1:  " << S.norm_operator_l1() << nl;
        LOG << "Norm Operator Max: " << S.norm_operator_max() << nl;
        LOG << nl;

        assert( S.norm_operator_l1()  == 11 );
        assert( S.norm_operator_max() == 10 );
        
        const auto GR = S.GerschgorinRow();
        const auto GC = S.GerschgorinColumn();
        LOG << "GerschgorinRow:    " << GR << nl;
        LOG << "GerschgorinColumn: " << GC << nl;

        assert( GR(0,0) ==  3 and GR(0,1) ==  6 );
        assert( GR(1,0) == -1 and GR(1,1) ==  3 );
        assert( GR(2,0) ==  1 and GR(2,1) ==  3 );
        assert( GR(3,0) ==  0 and GR(3,1) == 10 );

        assert( GC(0,0) ==  3 and GC(0,1) ==  4 );
        assert( GC(1,0) == -1 and GC(1,1) ==  5 );
        assert( GC(2,0) ==  1 and GC(2,1) == 10 );
        assert( GC(3,0) ==  0 and GC(3,1) ==  3 );
        
    }
    
    
    {
        LOG << "2. Transposing dense matrices" << nl;
        
        LOG << "2.A. Transpose of square matrices" << nl;
        
        for( int dim = 0; dim < 4; dim++ ) 
        {
            
            std::function<Float(int,int)> testmatrix = [](int r, int c) 
                                                        { return 2*r + c; };
            
            DenseMatrix A(dim, testmatrix);
            DenseMatrix B = A;
            
            auto At1 = Transpose( A );
            auto At2 = TransposeSquare( A );
            LOG << A << At1 << At2 << nl;

            for( int r = 0; r < dim; r++ )
            for( int c = 0; c < dim; c++ )
                assert( A(r,c) == At1(c,r) and A(r,c) == At2(c,r) );
            
            TransposeSquareInSitu( B );
            LOG << B << nl;
            
            for( int r = 0; r < dim; r++ )
            for( int c = 0; c < dim; c++ )
                assert( At1(r,c) == B(r,c) );
            
        }
        
        LOG << "2.B. Transpose of non-square matrices" << nl;
        
        for( int dimr = 0; dimr < 4; dimr++ ) 
        for( int dimc = 0; dimc < 4; dimc++ ) 
        {
            
            std::function<Float(int,int)> testmatrix = [](int r, int c) 
                                                        { return 2*r + c; };
            
            DenseMatrix A( dimr, dimc, testmatrix );
            const auto At = Transpose(A);
            
            for( int r = 0; r < dimr; r++ )
            for( int c = 0; c < dimc; c++ )
                assert( A(r,c) == At(c,r) );
            
            LOG << A << At << nl;
            //TransposeInSitu( A );
            // LOG << A << nl;

            // for( int r = 0; r < dimr; r++ )
            // for( int c = 0; c < dimc; c++ )
            //     assert( A(c,r) == At(c,r) );
            
        }
        
    }
    
    
    {
        
        LOG << "3. Cholesky decomposition of a fixed sample matrix" << nl;
        
        DenseMatrix A(3,3);
        
        A(0,0) =   4; A(0,1) =  12; A(0,2) = -16; 
        A(1,0) =  12; A(1,1) =  37; A(1,2) = -43; 
        A(2,0) = -16; A(2,1) = -43; A(2,2) =  98; 
        
        DenseMatrix L = CholeskyDecomposition( A );

        DenseMatrix ActualL( 3, 3 );
        ActualL(0,0) =  2; ActualL(0,1) = 0; ActualL(0,2) = 0; 
        ActualL(1,0) =  6; ActualL(1,1) = 1; ActualL(1,2) = 0; 
        ActualL(2,0) = -8; ActualL(2,1) = 5; ActualL(2,2) = 3; 
        
        //LOG << "Original matrix:\n" << A << nl;
        //LOG << "Factor matrix:\n" << L << nl;
        //LOG << "Product matrix:\n" << L * Transpose(L) << nl;

        assert( ( ActualL - L ).is_numerically_small() );
        assert( ( L * Transpose(L) - A ).is_numerically_small() );
    }
      
    {
    
        LOG << "4. Cholesky decompositions of a series of matrices" << nl;
        
        for( int dim = 0; dim <= 3; dim++ ) 
        {
            DenseMatrix A(dim);
            
            A.zero_matrix();
            for( int s = 0; s < dim; s++ )
            for( int t = 0; t < dim; t++ )
                A(s,t) = 3 * kronecker(s,t) - kronecker(s,t-1) - kronecker(s,t+1);
            
            DenseMatrix L = CholeskyDecomposition( A );
            
            LOG << "Original matrix: " << A << nl;
            
            LOG << "Factor matrix:   " << L << nl;
            
            LOG << "Product matrix:  " << L * Transpose(L) << nl;

            assert( ( L * Transpose(L) - A ).is_numerically_small() );
        }
        
    }


    
    {
        
        LOG << "5. Gauss Jordan algorithm with a fixed sample matrix" << nl;
    
        DenseMatrix A(4,4);
        
        A(0,0) = 1; A(0,1) = 2; A(0,2) = 3; A(0,3) = 5; 
        A(1,0) = 1; A(1,1) = 1; A(1,2) = 1; A(1,3) = 6; 
        A(2,0) = 3; A(2,1) = 3; A(2,2) = 1; A(2,3) = 2; 
        A(3,0) = 2; A(3,1) = 1; A(3,2) = 0; A(3,3) = 1; 
        
        LOG << "Determinant (default): " << Determinant(A) << nl;
        LOG << "Determinant (laplace): " << Determinant_laplaceexpansion(A) << nl;
        LOG << "Determinant (gauss):   " << Determinant_gauss(A) << nl;
      
        
//         A(0,0) = 4; A(0,1) = 1; A(0,2) = 0; 
//         A(1,0) = 0; A(1,1) = 1; A(1,2) = 0; 
//         A(2,0) = 1; A(2,1) = 0; A(2,2) = 1; 
        
        LOG << "Original matrix:\n" << A << nl;
        
        const DenseMatrix Ainv = GaussJordanInplace( A );
        
        LOG << "Proposed Inverse:\n" << Ainv << nl;
        
        LOG << "Product of both:\n" << Ainv * A << nl;

        assert( ( Ainv * A ).is_numerically_identity() );
        assert( ( A * Ainv ).is_numerically_identity() );
        
    }
    
    {
        LOG << "6. Gauss Jordan algorithm with several random integer matrices" << nl;
        for( int d = 0; d <= 8; d++ )
        for( int i = 0; i <  6; i++ )
        {
            
            DenseMatrix A(d);
            // LOG << A << nl;
            A.random_integer_matrix(-2,2);
            A += 30. * IdentityMatrix(d);
            LOG << A << nl;
            
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << "Determinant (default): " << det << nl;
            LOG << "Determinant (laplace): " << det_l << nl;
            LOG << "Determinant (gauss):   " << det_g << nl;
            LOG << "Determinant (bareiss): " << det_b << nl;
            
            // LOG << A << nl;

            Assert( is_numerically_one( det / det_l ), det, det_l );
            Assert( is_numerically_one( det / det_g ), det, det_g );
            Assert( is_numerically_one( det / det_b ), det, det_b );

            const auto Ainv = GaussJordanInplace(A);
            
            // LOG << Ainv * A << nl;

            assert( ( Ainv * A ).is_numerically_identity() );
            assert( ( A * Ainv ).is_numerically_identity() );
        }
    }

    {
        LOG << "7. Gauss Jordan algorithm with Hilbert matrices" << nl;
    
        for( int N = 1; N <= 6; N++ )
        {
            DenseMatrix C(N);
            for( int i = 0; i < N; i++ )
            for( int j = 0; j < N; j++ )
                C(i,j) = 1. / ( i+j+1 );
            
            auto Cinv = GaussJordanInplace(C);
            
            Assert( ( C * Cinv ).is_numerically_identity(), C * Cinv );
            Assert( ( Cinv * C ).is_numerically_identity(), Cinv * C );
        }
        
        
    }
    
    
    {
        LOG << "8. QR and QL Factorization" << nl;
    
        for( int dim = 0; dim <= 10; dim++ ){
            
            DenseMatrix A(dim,dim);
        
            A.zero_matrix();
            for( int s = 0; s < dim; s++ )
            for( int t = 0; t < dim; t++ )
                A(s,t) = 3 * kronecker(s,t) - kronecker(s,t-1) - kronecker(s,t+1);
                
            {
                DenseMatrix Q(dim,dim), R(dim,dim);
                
                QRFactorization( A, Q, R );
                
                DenseMatrix Rinv =   Inverse(R);
                DenseMatrix Qt   = Transpose(Q);
                
                assert( ( Q * R - A ).is_numerically_small() );
                assert( ( Q * Qt ).is_numerically_identity() );
                
                LOG << "Matrix A:\n" << A;
                LOG << "Matrix Q:\n" << Q;
                LOG << "Matrix R:\n" << R;
                LOG << "Matrix Q * R:\n" << Q * R;
                LOG << "Matrix Q^t:\n" << Qt;
                LOG << "Matrix Rinv:\n" << Rinv;
                LOG << "Matrix R * Rinv:\n" << R * Rinv;
                LOG << "Matrix Rinv * R:\n" << Rinv * R;
                LOG << "Matrix Q * Qinv:\n" << Q * Qt;
                LOG << "Matrix Qinv * Q:\n" << Qt * Q;
                LOG << "Matrix Rinv * Q^t * A:\n" << Rinv * Qt * A;
                LOG << "Matrix A * Rinv * Q^t:\n" << A * Rinv * Qt;
            }

            // if( /* DISABLES CODE */ (false) )
            {
                DenseMatrix Q(dim,dim), R(dim,dim);
                
                QRFactorization_via_Householder( A, Q, R );
                // TransposeSquareInSitu(Q);
                
                DenseMatrix Rinv =   Inverse(R);
                DenseMatrix Qt   = Transpose(Q);
                
                LOG << "Matrix A:\n" << A;
                LOG << "Matrix Q:\n" << Q;
                LOG << "Matrix R:\n" << R;
                LOG << "Matrix Q * R:\n" << Q * R;
                LOG << "Matrix Q^t:\n" << Qt;
                LOG << "Matrix Rinv:\n" << Rinv;
                LOG << "Matrix R * Rinv:\n" << R * Rinv;
                LOG << "Matrix Rinv * R:\n" << Rinv * R;
                LOG << "Matrix Q * Qinv:\n" << Q * Qt;
                LOG << "Matrix Qinv * Q:\n" << Qt * Q;
                LOG << "Matrix Rinv * Q^t * A:\n" << Rinv * Qt * A;
                LOG << "Matrix A * Rinv * Q^t:\n" << A * Rinv * Qt;
                
                assert( ( Q * R - A ).is_numerically_small() );
                assert( ( Q * Qt ).is_numerically_identity() );
            }

            {
                DenseMatrix L(dim,dim), Q(dim,dim);
                
                LQFactorization( A, L, Q );
                
                DenseMatrix Linv =   Inverse(L);
                DenseMatrix Qt   = Transpose(Q);

                assert( ( L * Q - A ).is_numerically_small() );
                assert( ( Q * Qt ).is_numerically_identity() );
                
                LOG << "Matrix A:\n" << A;
                LOG << "Matrix Q:\n" << Q;
                LOG << "Matrix L:\n" << L;
                LOG << "Matrix L * Q:\n" << L * Q;
                LOG << "Matrix Q^t:\n" << Qt;
                LOG << "Matrix Linv:\n" << Linv;
                LOG << "Matrix L * Linv:\n" << L * Linv;
                LOG << "Matrix Linv * L:\n" << Linv * L;
                LOG << "Matrix Q * Qinv:\n" << Q * Qt;
                LOG << "Matrix Qinv * Q:\n" << Qt * Q;
                LOG << "Matrix Q^t * Linv * A:\n" << Qt * Linv * A;
                LOG << "Matrix A * Qt * Linv:\n" << A * Qt * Linv;
            }

            {
                DenseMatrix Q1(dim,dim), R(dim,dim);
                QRFactorization( A, Q1, R );

                DenseMatrix L(dim,dim), Q2(dim,dim);
                LQFactorization( Transpose(A), L, Q2 );

                assert( ( Q1 * R - Transpose(Q2) * Transpose(L) ).is_numerically_small() );
                
            }


        }
        
        
    }
    
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
    
    
}
