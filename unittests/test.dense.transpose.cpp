

/**/

#include <iostream>
#include "../basic.hpp"
#include "../dense/dense.functions.hpp"


using namespace std;

int main()
{
    cout << "Unit test for transpoing dense matrices" << endl;
    
    std::cout.precision(5);
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    std::cout << std::showpos;
    
    std::function<Float(int,int)> testmatrix = [](int r, int c) 
                                                  { return 2*r + c; };
    
    std::function<Float(int,int)> hilbertmatrix = [](int r, int c) 
                                                  { return 1. / (r+c+1); };
    
    cout << "Transpose of square matrices" << endl;
    
    for( int dim = 0; dim < 4; dim++ ) 
    {
      
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
      
      DenseMatrix A( dimr, dimc, hilbertmatrix);
      
      cout << A << Transpose(A) << endl;
      TransposeInSitu( A );
      cout << A << endl;
      
    }
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
