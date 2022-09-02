

/**/

#include <ostream>
#include <iomanip>
#include "../../basic.hpp"
#include "../../dense/cholesky.hpp"
#include "../../dense/functions.hpp"
#include "../../dense/gaussjordan.hpp"
#include "../../dense/qr.factorization.hpp"
#include "../../dense/scalarfunctions.hpp"


using namespace std;

int main()
{
    LOG << "Unit Tests for Matrix Algorithms" << nl;
    
    // LOG << std::setprecision(5);
    // LOG << std::fixed << std::ios::floatfield; //setf( std::ios::fixed, std::ios::floatfield ); 
    // LOG << std::showpos;
    
    {
        LOG << "A. Unit test for scalar functions of matrices" << nl;
       
        DenseMatrix S( 4, 4 );
        S(0,0) =  3; S(0,1) =  0; S(0,2) = 6; S(0,3) =  0; 
        S(1,0) =  1; S(1,1) = -1; S(1,2) = 0; S(1,3) =  2; 
        S(2,0) = -1; S(2,1) =  1; S(2,2) = 1; S(2,3) = -1; 
        S(3,0) =  2; S(3,1) = -4; S(3,2) = 4; S(3,3) =  0; 
        
        
        LOG << S << nl;
        LOG << "Matrix trace:   " << MatrixTrace( S ) << nl;
        LOG << "Norm L1:        " << NormL1( S ) << nl;
        LOG << "Norm Frobenius: " << NormFrobenius( S ) << nl;
        LOG << "Norm Max:       " << NormMax( S ) << nl;
        
        Float p = 1.01;
        LOG << nl << "Norm Lp with p=" << p << ": " << NormLp( S, p ) << nl << nl;
        
        Float p1 = 100.0001;
        Float p2 = 1.00001;
        LOG << nl << "Row " << p1 << space
                     << "Col " << p2 << space
                     << NormRowCol( S, p1, p2 ) << nl;
        LOG << nl << "Col " << p1 << space
                     << "Row " << p2 << space 
                     << NormColRow( S, p1, p2 ) << nl;
        
        LOG << nl << "Row " << 1. << space
                     << "Col " << 1. << space
                     << NormRowCol( S, 1., 1. ) << nl;
        LOG << nl << "Col " << 1. << space
                     << "Row " << 1. << space 
                     << NormColRow( S, 1., 1. ) << nl;
        LOG << nl;
        
        LOG << nl << "Row " << 2. << space
                     << "Col " << 2. << space
                     << NormRowCol( S, 2., 2. ) << nl;
        LOG << nl << "Col " << 2. << space
                     << "Row " << 2. << space 
                     << NormColRow( S, 2., 2. ) << nl;
        LOG << nl;
        
        LOG << nl << "Row " << 20. << space
                     << "Col " << 20. << space
                     << NormRowCol( S, 20., 20. ) << nl;
        LOG << nl << "Col " << 20. << space
                     << "Row " << 20. << space 
                     << NormColRow( S, 20., 20. ) << nl;
        LOG << nl;
        
        
        
        LOG << "Norm Operator L1:  " << NormOperatorL1( S ) << nl;
        LOG << "Norm Operator Max: " << NormOperatorMax( S ) << nl;
        LOG << nl;
        
        LOG << "GerschgorinRow:    " << GerschgorinRow( S ) << nl;
        LOG << "GerschgorinColumn: " << GerschgorinColumn( S ) << nl;
        
    }
    
    
    {
        LOG << "B. Unit test for transpoing dense matrices" << nl;
        
        LOG << "Transpose of square matrices" << nl;
        
        for( int dim = 0; dim < 4; dim++ ) 
        {
            
            std::function<Float(int,int)> testmatrix = [](int r, int c) 
                                                        { return 2*r + c; };
            
            DenseMatrix A(dim, testmatrix);
            DenseMatrix B = A;
            
            LOG << A << Transpose( A ) << TransposeSquare( A ) << nl;
            TransposeInSitu( A );
            TransposeSquareInSitu( B );
            LOG << A << B << nl;
            
        }
        
        LOG << "Transpose of non-square matrices" << nl;
        
        for( int dimr = 0; dimr < 4; dimr++ ) 
        for( int dimc = 0; dimc < 4; dimc++ ) 
        {
            
            std::function<Float(int,int)> testmatrix = [](int r, int c) 
                                                        { return 2*r + c; };
            
            DenseMatrix A( dimr, dimc, testmatrix );
            
            LOG << A << Transpose(A) << nl;
            TransposeInSitu( A );
            LOG << A << nl;
            
        }
        
    }
    
    
    {
        LOG << "1. Unit test for matrix determinant, inverse, cofactor matrix" << nl;
        
        for( int dim = 1; dim < 6; dim++ ) 
        {
            DenseMatrix A = HilbertMatrix( dim );
            
            LOG << A << nl;
            
            LOG << "Determinant (default): " << Determinant(A) << nl;
            LOG << "Determinant (laplace): " << Determinant_laplaceexpansion(A) << nl;
            LOG << "Determinant (gauss):   " << Determinant_gauss(A) << nl;
            LOG << CofactorMatrix( A ) << nl;
            LOG << Inverse( A ) << nl;
            LOG << A * Inverse( A ) << nl;
            LOG << Inverse( A ) * A << nl;
            LOG << 1. / Determinant(A) - Determinant( Inverse(A) ) << nl;
            LOG << A - Inverse( Inverse(A) ) << nl;
        }
    }
    
    {
        
        LOG << "2. Unit Test for Cholesky decomposition" << nl;
        
        DenseMatrix A(3,3);
        
        A(0,0) =   4; A(0,1) =  12; A(0,2) = -16; 
        A(1,0) =  12; A(1,1) =  37; A(1,2) = -43; 
        A(2,0) = -16; A(2,1) = -43; A(2,2) =  98; 
        
        DenseMatrix L = CholeskyDecomposition( A );
        
        //LOG << "Original matrix:" << A << nl;
        //LOG << "Factor matrix:" << L << nl;
        //LOG << "Product matrix:" << L * Transpose(L) << nl;
        
    }
      
    {
    
        LOG << "3. Unit Test for Cholesky decomposition" << nl;
        
        int dim = 3;
        DenseMatrix A(dim);
        
        A.zeromatrix();
        for( int s = 0; s < dim; s++ )
        for( int t = 0; t < dim; t++ )
            A(s,t) = 3 * kronecker(s,t) - kronecker(s,t-1) - kronecker(s,t+1);
        
        DenseMatrix L = CholeskyDecomposition( A );
        
        //LOG << "Original matrix:" << A << nl;
        
        //LOG << "Factor matrix:" << L << nl;
        
        //LOG << "Product matrix:" << L * Transpose(L) << nl;
        
    }
    
    {
      
        LOG << "4. Unit test for matrix determinant, inverse, cofactor matrix" << nl;

        DenseMatrix A( 3 );
        A(0,0) = -2; A(0,1) =  2; A(0,2) = -3; 
        A(1,0) = -1; A(1,1) =  1; A(1,2) =  3; 
        A(2,0) =  2; A(2,1) =  0; A(2,2) = -1; 

        LOG << "Determinant (default): " << Determinant(A) << nl;
        LOG << "Determinant (laplace): " << Determinant_laplaceexpansion(A) << nl;
        LOG << "Determinant (gauss):   " << Determinant_gauss(A) << nl;
        LOG << CofactorMatrix( A ) << nl;
        LOG << Inverse( A ) << nl;
        LOG << A * Inverse( A ) << nl;
        LOG << Inverse( A ) * A << nl;
      
    }
    
    {
        
        LOG << "5. Unit Test for Gauss Jordan algorithm" << nl;
    
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
        
        LOG << "Original matrix:" << A << nl;
        
        DenseMatrix M = GaussJordanInplace( A );
//         DenseMatrix M = GaussJordanInplace( A );
        
        LOG << "Proposed Inverse:" << M << nl;
        
        LOG << "Product of both:" << M * A << nl;
        
    }
    
    {
        LOG << "6. Unit Test for Gauss Jordan algorithm" << nl;
        if(false)
        for( int i = 0; i < 6; i++ )
        {
            
            DenseMatrix A(8);
            A.randomintegermatrix(-2,2);
        
            LOG << "Determinant (default): " << Determinant(A) << nl;
            LOG << "Determinant (laplace): " << Determinant_laplaceexpansion(A) << nl;
            LOG << "Determinant (gauss):   " << Determinant_gauss(A) << nl;
        
            
            LOG << GaussJordanInplace(A) * A << nl;
        }
    }
    
    {
        LOG << "7. Unit Test for Gauss Jordan algorithm" << nl;
    
        int N = 14;
        DenseMatrix C(N);
        for( int i = 0; i < N; i++ )
        for( int j = 0; j < N; j++ )
            C(i,j) = 1. / ( i+j+1 );
        
        LOG << C * GaussJordanInplace(C,true) << nl;
        
        {
            auto Cinv = GaussJordanInplace(C);
            DenseMatrix I(N); I.unitmatrix();
            Float relaxation = Cinv.norm();
            
            for( int t = 0; t < 2000; t++ )
                Cinv = Cinv + 1.5/ ( std::log(N)) * ( I - C * Cinv );
            
            LOG << C * Cinv << nl;
        }
        
    }
    
    
    {
        LOG << "8. Unit Test for QR Factorization" << nl;
    
        const int dim = 4;
        DenseMatrix A(dim,dim);
    
        A.zeromatrix();
        for( int s = 0; s < dim; s++ )
        for( int t = 0; t < dim; t++ )
            A(s,t) = 3 * kronecker(s,t) - kronecker(s,t-1) - kronecker(s,t+1);
            
        {
            DenseMatrix Q(dim,dim), R(dim,dim);
            
            QRFactorization( A, Q, R );
            
            DenseMatrix Rinv =   Inverse(R);
            DenseMatrix Qt   = Transpose(Q);
            
            LOG << "Matrix A:" << A;
            LOG << "Matrix Q:" << Q;
            LOG << "Matrix R:" << R;
            LOG << "Matrix Q * R:" << Q * R;
            LOG << "Matrix Q^t:" << Qt;
            LOG << "Matrix Rinv:" << Rinv;
            LOG << "Matrix R * Rinv:" << R * Rinv;
            LOG << "Matrix Rinv * R:" << Rinv * R;
            LOG << "Matrix Q * Qinv:" << Q * Qt;
            LOG << "Matrix Qinv * Q:" << Qt * Q;
            LOG << "Matrix Rinv * Q^t * A:" << Rinv * Qt * A;
            LOG << "Matrix A * Rinv * Q^t:" << A * Rinv * Qt;
        }
        
    }
    
    {
        LOG << "9. Compare Determinats of Matrices with random coefficients" << nl;
        for( int t = 0; t < 7; t++ )
        for( int i = 0; i < 6; i++ )
        {
            
            DenseMatrix A(t);
            A.randomintegermatrix(10,80);
            
            LOG << A << Determinant_laplaceexpansion(A);
            
            LOG << "Determinant (laplace): " << Determinant_laplaceexpansion(A) << nl;
            LOG << "Determinant (gauss):   " << Determinant_gauss(A) << nl;
            LOG << "Determinant (default): " << Determinant(A) << nl;
            LOG << CofactorMatrix( A ) << nl;
            LOG << Inverse( A ) << nl;
            LOG << A * Inverse( A ) << nl;
            LOG << Inverse( A ) * A << nl;
        }
    }

   
    
    LOG << "Finished Unit Test" << nl;

    return 0;
    
    
}
