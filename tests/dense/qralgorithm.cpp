

/**/

#include <ostream>
#include "../../basic.hpp"
#include "../../dense/cholesky.hpp"
#include "../../dense/functions.hpp"
#include "../../dense/gaussjordan.hpp"
#include "../../dense/qr.factorization.hpp"
#include "../../dense/scalarfunctions.hpp"


using namespace std;

int main()
{
    LOG << "Unit Tests for QR Algorithm" << nl;
    
    {
        LOG << "8. Unit Test for QR Factorization" << nl;
    
        const int dim = 4;
        DenseMatrix A(dim,dim);

        A.zeromatrix();
        for( int s = 0; s < dim; s++ )
        for( int t = 0; t < dim; t++ )
            A(s,t) = 3 * kronecker(s,t) - kronecker(s,t-1) - kronecker(s,t+1);
            
        int repetitions = 4000;
        while ( repetitions --> 0 )
        {
            DenseMatrix Q(dim,dim), R(dim,dim);
            
            QRFactorization( A, Q, R );
            
            A = R * Q;
            
            if ( repetitions % 100 ) LOG << "Matrix A:" << A;
        }
        
    }
    
    
    
    LOG << "Finished Unit Test" << nl;

    return 0;
    
    
}
