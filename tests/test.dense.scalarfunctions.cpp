

/**/

#include <iostream>
#include "../basic.hpp"
#include "../dense/scalarfunctions.hpp"


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
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
