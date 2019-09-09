

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
        
        DenseMatrix A(3,3);
        
        A(0,0) = 1; A(0,1) = 2; A(0,2) = 3; 
        A(1,0) = 1; A(1,1) = 1; A(1,2) = 1; 
        A(2,0) = 3; A(2,1) = 3; A(2,2) = 1; 
        
//         A(0,0) = 4; A(0,1) = 1; A(0,2) = 0; 
//         A(1,0) = 0; A(1,1) = 1; A(1,2) = 0; 
//         A(2,0) = 1; A(2,1) = 0; A(2,2) = 1; 
        
        cout << "Original matrix:" << A << endl;
        
        DenseMatrix M = GaussJordan( A );
        
        cout << "Proposed Inverse:" << M << endl;
        
        cout << "Product of both:" << M * A << endl;
        
    }
      
    if(false){
    
        int dim = 3;
        DenseMatrix A(dim);
        
        A.zeromatrix();
        for( int s = 0; s < dim; s++ )
        for( int t = 0; t < dim; t++ )
            A(s,t) = 3 * kronecker(s,t) - kronecker(s,t-1) - kronecker(s,t+1);
        
        DenseMatrix M = GaussJordan( A );
        
        cout << "Original matrix:" << A << endl;
        
        cout << "Result:" << M << endl;
        
    }
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
