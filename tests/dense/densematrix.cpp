
#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Dense Matrix class" << nl;

    { 
        LOG << "Random matrix A of size 3x4, and 3*A" << nl;

        DenseMatrix A( 3, 4 );
        A.randommatrix();

        const auto Ax3 = 3 * A;

        LOG << A << nl;
        LOG << Ax3 << nl;

        for( int r = 0; r < 3; r++ )
        for( int c = 0; c < 4; c++ )
            assert( Ax3(r,c) == 3 * A(r,c) );

        LOG << "Random matrix B of size 3x4, and A+B" << nl;

        DenseMatrix B( 3, 4 );
        B.randommatrix();

        const auto C = A + B; 

        LOG << B << nl;
        LOG << C << nl;

        for( int r = 0; r < 3; r++ )
        for( int c = 0; c < 4; c++ )
            assert( C(r,c) == A(r,c) + B(r,c) );

    }

    {
        LOG << "Unit Matrices of size 3x3 and 4x4" << nl;

        DenseMatrix A( 3, 4 );
        A.randommatrix();

        DenseMatrix I3(3,3);
        I3.unitmatrix();
        DenseMatrix I4(4,4);
        I4.unitmatrix();
        LOG << I3 << I4 << nl;
        
        LOG << "I3 * A and A * I4" << nl;

        const auto I3xA = I3 * A;
        const auto AxI4 = A * I4;

        LOG << I3xA << nl;
        LOG << AxI4 << nl;

        for( int r = 0; r < 3; r++ )
        for( int c = 0; c < 4; c++ )
            assert( A(r,c) == I3xA(r,c) and A(r,c) == AxI4(r,c) );
        
        LOG << "5 * I3 and (5*I3)* A" << nl;

        const auto S5 = 5. * I3;
        const auto S5xA = operator*( S5, A );
        LOG << S5 << nl;
        LOG << S5xA << nl;

        for( int r = 0; r < 3; r++ )
        for( int c = 0; c < 4; c++ )
            assert( S5xA(r,c) == 5. * A(r,c) );

    }
    

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
