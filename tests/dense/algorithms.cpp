

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../dense/cholesky.hpp"
#include "../../dense/functions.hpp"
#include "../../dense/gaussjordan.hpp"
#include "../../dense/qr.factorization.hpp"
#include "../../dense/scalarfunctions.hpp"


using namespace std;

int main()
{
    cout << "Unit Tests for Matrix Algorithms" << endl;
    
    std::cout.precision(5);
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    std::cout << std::showpos;
    
    {
        
        cout << "Unit Test for Cholesky decomposition" << endl;
        
        DenseMatrix A(3,3);
        
        A(0,0) =   4; A(0,1) =  12; A(0,2) = -16; 
        A(1,0) =  12; A(1,1) =  37; A(1,2) = -43; 
        A(2,0) = -16; A(2,1) = -43; A(2,2) =  98; 
        
        DenseMatrix L = CholeskyDecomposition( A );
        
        cout << "Original matrix:" << A << endl;
        
        cout << "Factor matrix:" << L << endl;
        
        cout << "Product matrix:" << L * Transpose(L) << endl;
        
    }
      
    {
    
        cout << "Unit Test for Cholesky decomposition" << endl;
        
        int dim = 3;
        DenseMatrix A(dim);
        
        A.zeromatrix();
        for( int s = 0; s < dim; s++ )
        for( int t = 0; t < dim; t++ )
            A(s,t) = 3 * kronecker(s,t) - kronecker(s,t-1) - kronecker(s,t+1);
        
        DenseMatrix L = CholeskyDecomposition( A );
        
        cout << "Original matrix:" << A << endl;
        
        cout << "Factor matrix:" << L << endl;
        
        cout << "Product matrix:" << L * Transpose(L) << endl;
        
    }
    
    for( int dim = 0; dim < 6; dim++ ) 
    {
      
        cout << "Unit test for matrix determinant, inverse, cofactor matrix" << endl;

        DenseMatrix A = HilbertMatrix( dim );
        
        cout << A << endl;
        
        cout << "Determinant (default): " << Determinant(A) << endl;
        cout << "Determinant (laplace): " << Determinant_laplaceexpansion(A) << endl;
        cout << "Determinant (gauss):   " << Determinant_gauss(A) << endl;
        cout << CofactorMatrix( A ) << endl;
        cout << Inverse( A ) << endl;
        cout << A * Inverse( A ) << endl;
        cout << Inverse( A ) * A << endl;
        cout << 1. / Determinant(A) - Determinant( Inverse(A) ) << endl;
        cout << A - Inverse( Inverse(A) ) << endl;
            
    }
    
    {
      
        cout << "Unit test for matrix determinant, inverse, cofactor matrix" << endl;

        DenseMatrix A( 3 );
        A(0,0) = -2; A(0,1) =  2; A(0,2) = -3; 
        A(1,0) = -1; A(1,1) =  1; A(1,2) =  3; 
        A(2,0) =  2; A(2,1) =  0; A(2,2) = -1; 

        cout << "Determinant (default): " << Determinant(A) << endl;
        cout << "Determinant (laplace): " << Determinant_laplaceexpansion(A) << endl;
        cout << "Determinant (gauss):   " << Determinant_gauss(A) << endl;
        cout << CofactorMatrix( A ) << endl;
        cout << Inverse( A ) << endl;
        cout << A * Inverse( A ) << endl;
        cout << Inverse( A ) * A << endl;
      
    }
    
    cout << "Compare Determinats of Matrices with random coefficients" << endl;
    for( int t = 0; t < 7; t++ )
    for( int i = 0; i < 6; i++ )
    {
        
        DenseMatrix A(t);
        A.randomintegermatrix(-5,5);
        
        cout << "Determinant " << t << " (default/gauss/laplace): "
             << Determinant(A) << space 
             << Determinant_gauss(A) << space 
             << Determinant_laplaceexpansion(A) << space 
             << 1. / Determinant(A) - Determinant( Inverse(A) ) << space
             << endl;
            
        if( absolute( Determinant_gauss(A) - Determinant_laplaceexpansion(A) ) > 0.000001
            or
            absolute( Determinant_gauss(A) - Determinant(A) ) > 0.000001
        ) {
            cout << A; return 1;
        }
        
    }

    {
        
        cout << "Unit Test for Gauss Jordan algorithm" << endl;
    
        DenseMatrix A(4,4);
        
        A(0,0) = 1; A(0,1) = 2; A(0,2) = 3; A(0,3) = 5; 
        A(1,0) = 1; A(1,1) = 1; A(1,2) = 1; A(1,3) = 6; 
        A(2,0) = 3; A(2,1) = 3; A(2,2) = 1; A(2,3) = 2; 
        A(3,0) = 2; A(3,1) = 1; A(3,2) = 0; A(3,3) = 1; 
        
        cout << "Determinant (default): " << Determinant(A) << endl;
        cout << "Determinant (laplace): " << Determinant_laplaceexpansion(A) << endl;
        cout << "Determinant (gauss):   " << Determinant_gauss(A) << endl;
      
        
//         A(0,0) = 4; A(0,1) = 1; A(0,2) = 0; 
//         A(1,0) = 0; A(1,1) = 1; A(1,2) = 0; 
//         A(2,0) = 1; A(2,1) = 0; A(2,2) = 1; 
        
        cout << "Original matrix:" << A << endl;
        
        DenseMatrix M = GaussJordanInplace( A );
//         DenseMatrix M = GaussJordanInplace( A );
        
        cout << "Proposed Inverse:" << M << endl;
        
        cout << "Product of both:" << M * A << endl;
        
    }
      
    for( int i = 0; i < 6; i++ )
    {
        cout << "Unit Test for Gauss Jordan algorithm" << endl;
    
        DenseMatrix A(8);
        A.randomintegermatrix(-2,2);
      
        cout << "Determinant (default): " << Determinant(A) << endl;
        cout << "Determinant (laplace): " << Determinant_laplaceexpansion(A) << endl;
        cout << "Determinant (gauss):   " << Determinant_gauss(A) << endl;
      
        
        cout << GaussJordanInplace(A) * A << nl;
    }

    
    {
        
        cout << "Unit Test for Gauss Jordan algorithm" << endl;
    
        int N = 14;
        DenseMatrix C(N);
        for( int i = 0; i < N; i++ )
        for( int j = 0; j < N; j++ )
            C(i,j) = 1. / ( i+j+1 );
        
        cout << C * GaussJordanInplace(C,true) << nl;
        
        {
            auto Cinv = GaussJordanInplace(C);
            DenseMatrix I(N); I.unitmatrix();
            Float relaxation = Cinv.norm();
            
            for( int t = 0; t < 2000; t++ )
                Cinv = Cinv + 1.5/ ( std::log(N)) * ( I - C * Cinv );
            
            cout << C * Cinv << nl;
        }
        
    }
    
    
    {
        
        cout << "Unit Test for QR Factorization" << endl;
    
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
            
            cout << "Matrix A:" << A;
            cout << "Matrix Q:" << Q;
            cout << "Matrix R:" << R;
            cout << "Matrix Q * R:" << Q * R;
            cout << "Matrix Q^t:" << Qt;
            cout << "Matrix Rinv:" << Rinv;
            cout << "Matrix R * Rinv:" << R * Rinv;
            cout << "Matrix Rinv * R:" << Rinv * R;
            cout << "Matrix Q * Qinv:" << Q * Qt;
            cout << "Matrix Qinv * Q:" << Qt * Q;
            cout << "Matrix Rinv * Q^t * A:" << Rinv * Qt * A;
            cout << "Matrix A * Rinv * Q^t:" << A * Rinv * Qt;
        }
        
    }
    
    
    
    
    {
        cout << "Unit test for scalar functions of matrices" << endl;
       
        DenseMatrix S( 4, 4 );
        S(0,0) =  3; S(0,1) =  0; S(0,2) = 6; S(0,3) =  0; 
        S(1,0) =  1; S(1,1) = -1; S(1,2) = 0; S(1,3) =  2; 
        S(2,0) = -1; S(2,1) =  1; S(2,2) = 1; S(2,3) = -1; 
        S(3,0) =  2; S(3,1) = -4; S(3,2) = 4; S(3,3) =  0; 
        
        
        cout << S << endl;
        cout << "Matrix trace:   " << MatrixTrace( S ) << endl;
        cout << "Norm L1:        " << NormL1( S ) << endl;
        cout << "Norm Frobenius: " << NormFrobenius( S ) << endl;
        cout << "Norm Max:       " << NormMax( S ) << endl;
        
        Float p = 1.01;
        cout << endl << "Norm Lp with p=" << p << ": " << NormLp( S, p ) << endl << endl;
        
        Float p1 = 100.0001;
        Float p2 = 1.00001;
        cout << endl << "Row " << p1 << space
                     << "Col " << p2 << space
                     << NormRowCol( S, p1, p2 ) << endl;
        cout << endl << "Col " << p1 << space
                     << "Row " << p2 << space 
                     << NormColRow( S, p1, p2 ) << endl;
        cout << endl;
        
        cout << endl << "Row " << 1. << space
                     << "Col " << 1. << space
                     << NormRowCol( S, 1., 1. ) << endl;
        cout << endl << "Col " << 1. << space
                     << "Row " << 1. << space 
                     << NormColRow( S, 1., 1. ) << endl;
        cout << endl;
        
        cout << endl << "Row " << 2. << space
                     << "Col " << 2. << space
                     << NormRowCol( S, 2., 2. ) << endl;
        cout << endl << "Col " << 2. << space
                     << "Row " << 2. << space 
                     << NormColRow( S, 2., 2. ) << endl;
        cout << endl;
        
        cout << endl << "Row " << 20. << space
                     << "Col " << 20. << space
                     << NormRowCol( S, 20., 20. ) << endl;
        cout << endl << "Col " << 20. << space
                     << "Row " << 20. << space 
                     << NormColRow( S, 20., 20. ) << endl;
        cout << endl;
        
        
        
        cout << "Norm Operator L1:  " << NormOperatorL1( S ) << endl;
        cout << "Norm Operator Max: " << NormOperatorMax( S ) << endl;
        cout << endl;
        
        cout << "GerschgorinRow:    " << GerschgorinRow( S ) << endl;
        cout << "GerschgorinColumn: " << GerschgorinColumn( S ) << endl;
        
    }
    
    
    {
        cout << "Unit test for transpoing dense matrices" << endl;
        
        cout << "Transpose of square matrices" << endl;
        
        for( int dim = 0; dim < 4; dim++ ) 
        {
            
            std::function<Float(int,int)> testmatrix = [](int r, int c) 
                                                        { return 2*r + c; };
            
            DenseMatrix A(dim, testmatrix);
            DenseMatrix B = A;
            
            cout << A << Transpose( A ) << TransposeSquare( A ) << endl;
            TransposeInSitu( A );
            TransposeSquareInSitu( B );
            cout << A << B << endl;
            
        }
        
        cout << "Transpose of non-square matrices" << endl;
        
        for( int dimr = 0; dimr < 4; dimr++ ) 
        for( int dimc = 0; dimc < 4; dimc++ ) 
        {
            
            std::function<Float(int,int)> testmatrix = [](int r, int c) 
                                                        { return 2*r + c; };
            
            DenseMatrix A( dimr, dimc, testmatrix );
            
            cout << A << Transpose(A) << endl;
            TransposeInSitu( A );
            cout << A << endl;
            
        }
        
    }




    
    
    cout << "Finished Unit Test" << endl;

    return 0;
    
    
}
