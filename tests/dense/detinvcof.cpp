

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
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
