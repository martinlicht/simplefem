
#include "../../basic.hpp"
#include "../../dense/cholesky.hpp"
#include "../../dense/functions.hpp"
#include "../../dense/gaussjordan.hpp"
#include "../../dense/qr.factorization.hpp"
// #include "../../dense/scalarfunctions.hpp"


int main( int argc, char *argv[] )
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
            
            if( repetitions % 1000 == 0 ) LOG << "Matrix A:" << A;
        }
        
    }
    
    
    
    {
        LOG << "Solving Least-Squares problem with QR factorization" << nl;
    
        const int cols =  4;
        const int rows = 40;
        
        DenseMatrix A(rows,cols);
        A.randommatrix();
        
        FloatVector b( A.getdimout() );
        b.random();
        
        auto x = SolveOverconstrained( A, b );
        
        auto r = b - A * x;
        
        LOG << "residual: " << r.norm() << nl;
        LOG << "orthogonality: " << (Transpose(A) * r).norm() << nl;
        
    }
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
    
    
}
