

/**/

#include <iostream>
#include "../basic.hpp"
#include "../dense/functions.hpp"
#include "../dense/scalarfunctions.hpp"
#include "../dense/iterativeinverse.hpp"


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
    
//     cout << "Produce Hilbert matrices" << endl;
    
    for( int dim = 0; false && dim < 6; dim++ ) 
    {
      
      DenseMatrix A( dim, hilbertmatrix );
      DenseMatrix X( dim );
      
      cout << A << endl;
      
      cout << Inverse( A ) << endl;
      newtoniteration( X, A, 10 );
      cout << Inverse( A ) << endl;
      cout << X << endl;
      
    }
    
    {
      
      DenseMatrix A( 3 );
      A(0,0) = -2; A(0,1) =  2; A(0,2) = -3; 
      A(1,0) = -1; A(1,1) =  1; A(1,2) =  3; 
      A(2,0) =  2; A(2,1) =  0; A(2,2) = -1; 
      
      DenseMatrix X( 3 );
      
      Float alpha = EigenvalueEstimate( A );
  
      DenseMatrix I(3);
      I.unitmatrix();
      
//       X = ( 1. / ( alpha * alpha ) ) * A;
//       X = 0.0000001 * A;
      
      cout << alpha << endl;
      
//       for( int i = 0; i < 25; i++ )
//         X = ( 2 * X ) - ( X * ( A * X ) );
      
      for( int i = 0; i < 100000000; i++ )
        X = X - 0.000001 * ( A * X - I ) ;
      
      cout << X * A << endl;
      cout << A * Inverse( A ) << endl;
      
    }
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
