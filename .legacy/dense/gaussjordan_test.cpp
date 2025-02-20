

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../dense/functions.hpp"
#include "../../dense/gaussjordan.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for Gauss Jordan algorithm" << endl;
    
    std::cout.precision(5);
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    std::cout << std::showpos;
    
    {
        
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
        DenseMatrix A(8);
        A.randomintegermatrix(-2,2);
      
        cout << "Determinant (default): " << Determinant(A) << endl;
        cout << "Determinant (laplace): " << Determinant_laplaceexpansion(A) << endl;
        cout << "Determinant (gauss):   " << Determinant_gauss(A) << endl;
      
        
        cout << GaussJordanInplace(A) * A << nl;
    }

    
    {
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
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
