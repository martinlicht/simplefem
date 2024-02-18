

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../dense/functions.hpp"
#include "../../dense/simplesolver.hpp"
#include "../../dense/qr.factorization.hpp"
#include "../../dense/gaussjordan.hpp"
#include "../../dense/cholesky.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for Dense Matrix Module" << endl;
    
    std::cout.precision(5);
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    std::cout << std::showpos;
    
    {
        
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
    
    cout << "Finished Unit Test" << endl;
    
    return 0;
}
