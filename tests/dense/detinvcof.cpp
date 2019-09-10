

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../dense/functions.hpp"


using namespace std;

int main()
{
    cout << "Unit test for matrix determinant, inverse, cofactor matrix" << endl;
    
    std::cout.precision(10);
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    std::cout << std::showpos;
    
    std::function<Float(int,int)> testmatrix = [](int r, int c) 
                                                  { return 2*r + c; };
    
    std::function<Float(int,int)> hilbertmatrix = [](int r, int c) 
                                                  { return 1. / (r+c+1); };
    
    cout << "Produce Hilbert matrices" << endl;
    
    for( int dim = 0; dim < 6; dim++ ) 
    {
      
      DenseMatrix A( dim, hilbertmatrix );
      
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
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
