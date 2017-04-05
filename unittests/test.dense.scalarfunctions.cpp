

/**/

#include <iostream>
#include "../basic.hpp"
#include "../dense/dense.scalarfunctions.hpp"


using namespace std;

int main()
{
    cout << "Unit test for scalar functions of matrices" << endl;
    
    std::cout.precision(10);
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    std::cout << std::showpos;
    
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
    cout << endl << "Norm Lp with p=" << p << ": " << NormLp( S, 1.01 ) << endl << endl;
    
    cout << "Norm Operator L1:  " << NormOperatorL1( S ) << endl;
    cout << "Norm Operator Max: " << NormOperatorMax( S ) << endl;
    
    cout << "GerschgorinRow:    " << GerschgorinRow( S ) << endl;
    cout << "GerschgorinColumn: " << GerschgorinColumn( S ) << endl;
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
