
#include <limits>

#include "../../base/include.hpp"
#include "../../utility/random.hpp"
#include "../../dense/densematrix.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test: Dense Matrix class" << nl;

    { 
        LOG << "Random matrix A of size 3x4, and 3*A" << nl;

        DenseMatrix A( 3, 4 );
        A.random_matrix();

        const auto Ax3 = 3 * A;

        LOG << A << nl;
        LOG << Ax3 << nl;

        for( int r = 0; r < 3; r++ )
        for( int c = 0; c < 4; c++ )
            assert( Ax3(r,c) == 3 * A(r,c) );

        LOG << "Random matrix B of size 3x4, and A+B" << nl;

        DenseMatrix B( 3, 4 );
        B.random_matrix();

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
        A.random_matrix();

        DenseMatrix I3(3,3);
        I3.identity_matrix();
        DenseMatrix I4(4,4);
        I4.identity_matrix();
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


    {
        LOG << "Matrix properties" << nl;
    
        const int max_m = 6;
        const int max_n = 6;
        for( int m = 0; m <= max_m; m++ ) 
        for( int n = 0; n <= max_n; n++ ) 
        {
            // Create an m x n DenseMatrix
            DenseMatrix mat(m, n, 0.);

            LOG << m << space << n << nl;

            assert( mat.is_square() == (m == n) );
            
            assert( mat.is_diagonal() == (m == n) );
            assert( mat.is_antisymmetric() == (m == n) );
            assert( mat.is_symmetric() == (m == n) );
            
            assert( mat.is_lower_left_triangular() == (m == n) );
            assert( mat.is_upper_left_triangular() == (m == n) );
            assert( mat.is_lower_right_triangular() == (m == n) );
            assert( mat.is_upper_right_triangular() == (m == n) );

            assert( mat.is_finite() );
            
            if( mat.getdimout() == 0 or mat.getdimin() == 0 )
            {
                assert( mat.is_positive() );
                assert( mat.is_negative() );
                assert( mat.is_nonpositive() );
                assert( mat.is_nonnegative() );
            }

            if( mat.getdimout() == 0 or mat.getdimin() == 0 ) continue;
            

            // Test is_zero()
            assert( mat.is_zero() == true );  // Empty DenseMatrix is treated as zero

            mat( random_integer() % mat.getdimout(), random_integer() % mat.getdimin() ) = Constants::feigenbaum_first;
            assert( mat.is_zero() == false );  // Empty DenseMatrix is treated as zero


            // Test is_finite()
            assert( mat.is_finite() == true );
            for( int r = 0; r < m; r++ ) 
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = (r==c)?( std::numeric_limits<double>::infinity()):notanumber;
            }
            assert( mat.is_finite() == false );



            // Test positivity and negativity
            for( int r = 0; r < m; r++ ) 
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = ( 1);
            }
            assert( mat.is_positive() == true );
            assert( mat.is_nonnegative() == true );
            assert( mat.is_negative() == false );
            assert( mat.is_nonpositive() == false );

            for( int r = 0; r < m; r++ ) 
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = (((r+c)%2)?1:0);
            }
            assert( mat.is_positive() == false );
            assert( mat.is_nonnegative() == true );
            assert( mat.is_negative() == false );
            assert( mat.is_nonpositive() == (m==1 and n==1));

            for( int r = 0; r < m; r++ ) 
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = ( -1);
            }
            assert( mat.is_positive() == false );
            assert( mat.is_nonnegative() == false );
            assert( mat.is_negative() == true );
            assert( mat.is_nonpositive() == true );

            for( int r = 0; r < m; r++ ) 
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = (((r+c)%2)?-1:0);
            }
            assert( mat.is_positive() == false );
            assert( mat.is_nonnegative() == (m==1 and n==1));
            assert( mat.is_negative() == false );
            assert( mat.is_nonpositive() == true );




            // Test is_symmetric() and is_antisymmetric()
            for( int r = 0; r < m; r++ ) 
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = r+c;
            }
            assert( mat.is_symmetric()     == (m==n) );
            assert( mat.is_antisymmetric() == (m==1 and n==1) );
            
            for( int r = 0; r < m; r++ ) 
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = (r+c) * (r!=c?(r>c?-1:1):0);
            }
            assert( mat.is_symmetric()     == (m==1 and n==1) );
            assert( mat.is_antisymmetric() == (m==n) );
            
            
            
            // Test is_diagonal()
            for( int r = 0; r < m; r++ ) 
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = ( (r == c) ? 1 : 0);
            }
            assert( mat.is_diagonal()               == (m == n)        );
            assert( mat.is_upper_left_triangular()  == (m==1 and n==1) );
            assert( mat.is_upper_right_triangular() == (m == n)        );
            assert( mat.is_lower_left_triangular()  == (m == n)        );
            assert( mat.is_lower_right_triangular() == (m==1 and n==1) );

            // Test upper left triangular properties
            for( int r = 0; r < m; r++ ) 
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = ( (r <= mat.getdimin()-1-c) ? 1 : 0);  // upper left
            }
            assert( mat.is_diagonal()               == (m==1 and n==1) );
            assert( mat.is_upper_left_triangular()  == (m == n)        );
            assert( mat.is_upper_right_triangular() == (m==1 and n==1) );
            assert( mat.is_lower_left_triangular()  == (m==1 and n==1) );
            assert( mat.is_lower_right_triangular() == (m==1 and n==1) );

            // Test upper right triangular properties
            for( int r = 0; r < m; r++ )
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = ( (r <= c) ? 2 : 0);  // upper right 
            }
            assert( mat.is_diagonal()               == (m==1 and n==1) );
            assert( mat.is_upper_left_triangular()  == (m==1 and n==1) );
            assert( mat.is_upper_right_triangular() == (m == n)        );
            assert( mat.is_lower_left_triangular()  == (m==1 and n==1) );
            assert( mat.is_lower_right_triangular() == (m==1 and n==1) );

            // Test lower left triangular properties
            for( int r = 0; r < m; r++ ) 
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = ( (r >= c) ? 3 : 0);  // lower left 
            }
            assert( mat.is_diagonal()               == (m==1 and n==1) );
            assert( mat.is_upper_left_triangular()  == (m==1 and n==1) );
            assert( mat.is_upper_right_triangular() == (m==1 and n==1) );
            assert( mat.is_lower_left_triangular()  == (m == n)        );
            assert( mat.is_lower_right_triangular() == (m==1 and n==1) );

            // Test lower right triangular properties
            for( int r = 0; r < m; r++ )
            for( int c = 0; c < n; c++ ) 
            {
                mat( r, c ) = ( (r >= mat.getdimin()-1-c) ? 4 : 0);  // lower right 
            }
            assert( mat.is_diagonal()               == (m==1 and n==1) );
            assert( mat.is_upper_left_triangular()  == (m==1 and n==1) );
            assert( mat.is_upper_right_triangular() == (m==1 and n==1) );
            assert( mat.is_lower_left_triangular()  == (m==1 and n==1) );
            assert( mat.is_lower_right_triangular() == (m==n)          );

        }


    }
    

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
