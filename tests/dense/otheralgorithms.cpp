
#include "../../basic.hpp"
#include "../../dense/factorization.hpp"
#include "../../dense/functions.hpp"
#include "../../dense/simplesolver.hpp"
// #include "../../dense/scalarfunctions.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Tests for Matrix Algorithms: inverses, determinants, cofactors" << nl;
    
    {
      
        LOG << "1. Matrix determinant, inverse, cofactor matrix of a fixed sample matrix" << nl;

        DenseMatrix A( 3 );
        A(0,0) = -2; A(0,1) =  2; A(0,2) = -3; 
        A(1,0) = -1; A(1,1) =  1; A(1,2) =  3; 
        A(2,0) =  2; A(2,1) =  0; A(2,2) = -1; 

        Float det = Determinant(A);
        Float det_l = Determinant_laplaceexpansion(A);
        Float det_g = Determinant_gauss(A);
        Float det_b = Determinant_bareiss(A);

        LOG << "Determinant (default): " << det << nl;
        LOG << "Determinant (laplace): " << det_l << nl;
        LOG << "Determinant (gauss):   " << det_g << nl;
        LOG << "Determinant (bareiss): " << det_b << nl;

        assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );

        const auto Acof = CofactorMatrix( A );
        const auto Ainv = Inverse( A );

        LOG << Acof << nl;
        LOG << Ainv << nl;
        LOG << A * Ainv << nl;
        LOG << Ainv * A << nl;

        assert( ( Ainv * A ).is_numerically_identity() );
        assert( ( A * Ainv ).is_numerically_identity() );
      
    }
    
    {
        LOG << "2. Matrix determinant, inverse, cofactor matrix with Hilbert matrices" << nl;
        
        for( int dim = 0; dim < 6; dim++ ) 
        {
            DenseMatrix A = HilbertMatrix( dim );
            
            LOG << A << nl;
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << "Determinant (default): " << det << nl;
            LOG << "Determinant (laplace): " << det_l << nl;
            LOG << "Determinant (gauss):   " << det_g << nl;
            LOG << "Determinant (bareiss): " << det_b << nl;

            assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );

            const auto Acof = CofactorMatrix( A );
            const auto Ainv = Inverse( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;
            LOG << 1. / Determinant(A) - Determinant( Ainv ) << nl;
            LOG << A - Inverse( Ainv ) << nl;

            Assert( ( Ainv * A ).is_numerically_identity(), Ainv, A, Ainv*A );
            Assert( ( A * Ainv ).is_numerically_identity(), A, Ainv, Ainv*A );
        }
    }
    
    
    {
        LOG << "3. Compare Determinants of Matrices with random coefficients" << nl;
        for( int t = 0; t < 7; t++ )
        for( int i = 0; i < 6; i++ )
        {
            
            DenseMatrix A(t);
            A.random_integer_matrix(10,80);
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << "Determinant (default): " << det << nl;
            LOG << "Determinant (laplace): " << det_l << nl;
            LOG << "Determinant (gauss):   " << det_g << nl;
            LOG << "Determinant (bareiss): " << det_b << nl;

            // LOGPRINTF( "%.17e %.17e %.17e %.17e \n", det, det_l, det_g, det_b );

            assert( is_numerically_one( det / det_l ) && is_numerically_one( det / det_g ) && is_numerically_one( det / det_b ) );
            
            const auto Acof = CofactorMatrix( A );
            const auto Ainv = Inverse( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;

            assert( ( Ainv * A ).is_numerically_identity() );
            assert( ( A * Ainv ).is_numerically_identity() );

            // LOG << "with Newton-Schulz \n";
            // auto Ainvp = Ainv; NewtonSchulz( A, Ainvp, 20 );

            // LOG << A * Ainvp << nl;
            // LOG << Ainvp * A << nl;

            // auto R1 = ( A * Ainvp - IdentityMatrix(t) ).norm() / ( A * Ainv - IdentityMatrix(t) ).norm();
            // auto R2 = ( Ainvp * A - IdentityMatrix(t) ).norm() / ( Ainv * A - IdentityMatrix(t) ).norm(); 
            
            // LOG << "R1 = " << R1 << space << "R2 = " << R2 << nl;
            // if( t >= 5 ) assert( R1 < 1.0 and R2 < 1.0 );

            // assert( ( Ainvp * A ).is_numerically_identity() );
            // assert( ( A * Ainvp ).is_numerically_identity() );
        }
    }
    
    {
        LOG << "4. Inverse and determinants of unit matrices" << nl;
        for( int t = 0; t < 7; t++ )
        {
            
            DenseMatrix A = 3. * IdentityMatrix(t);
            
            // LOG << A << nl;
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << "Determinant (default): " << det << nl;
            LOG << "Determinant (laplace): " << det_l << nl;
            LOG << "Determinant (gauss):   " << det_g << nl;
            LOG << "Determinant (bareiss): " << det_b << nl;

            assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );
            
            const auto Acof = CofactorMatrix( A );
            const auto Ainv = GaussJordanInplace( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;

            assert( ( Ainv * A ).is_numerically_identity() );
            assert( ( A * Ainv ).is_numerically_identity() );
        }
    }

    
    // if(false)
    {
        LOG << "5. Compare Determinants of random orthogonal Matrices" << nl;
        for( int t = 0; t < 7; t++ )
        for( int i = 0; i < 6; i++ )
        {
            
            DenseMatrix A(t);
            A.random_orthogonal_matrix();

            LOG << A << nl;
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << "Determinant (default): " << det << nl;
            LOG << "Determinant (laplace): " << det_l << nl;
            LOG << "Determinant (gauss):   " << det_g << nl;
            LOG << "Determinant (bareiss): " << det_b << nl;

            assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );
            
            {
                const int dim = t;
                DenseMatrix Q(dim,dim), R(dim,dim);
                QRFactorization( A, Q, R );
                Float detR = UpperTriangularDeterminant(R); Float detQ = Determinant(Q);
                LOG << "Determinant (QR): " << detR << nl;
                Assert( is_numerically_one( detQ ) or is_numerically_one( -detQ ) );
                Assert( is_numerically_close( det, detQ * detR ) );

            }
            
            const auto Acof = CofactorMatrix( A );
            const auto Ainv = Inverse( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;
            
            assert( ( Ainv * A ).is_numerically_identity() );
            assert( ( A * Ainv ).is_numerically_identity() );

            // LOG << "with Newton-Schulz \n";
            // auto Ainvp = Ainv; NewtonSchulz( A, Ainvp, 200 );

            // LOG << A * Ainvp << nl;
            // LOG << Ainvp * A << nl;

            // auto R1 = ( A * Ainvp - IdentityMatrix(t) ).norm() / ( A * Ainv - IdentityMatrix(t) ).norm();
            // auto R2 = ( Ainvp * A - IdentityMatrix(t) ).norm() / ( Ainv * A - IdentityMatrix(t) ).norm(); 
            
            // LOG << "R1 = " << R1 << space << "R2 = " << R2 << nl;
            // if( t >= 5 ) assert( R1 < 1.0 and R2 < 1.0 );

            // assert( ( Ainvp * A ).is_numerically_identity() );
            // assert( ( A * Ainvp ).is_numerically_identity() );
        }
    }
    
    
    {
        LOG << "6. Compare Determinants of random ill-conditioned 2x2 matrices" << nl;
        for( int t = 3; t < 5; t++ )
        {
            
            DenseMatrix A(2,2);
            A(0,0) = power_numerical( 10, t );
            A(0,1) = power_numerical( 10, t ) - 1;
            A(1,0) = power_numerical( 10, t ) - 1;
            A(1,1) = power_numerical( 10, t ) - 2;
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << t << " Determinant (default): " << det << nl;
            LOG << t << " Determinant (laplace): " << det_l << nl;
            LOG << t << " Determinant (gauss):   " << det_g << nl;
            LOG << t << " Determinant (bareiss): " << det_b << nl;

            assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );
            
            {
                const int dim = 2;
                DenseMatrix Q(dim,dim), R(dim,dim);
                QRFactorization( A, Q, R );
                Float detR = UpperTriangularDeterminant(R); Float detQ = Determinant(Q);
                LOG << "Determinant (QR): " << detR << nl;
                Assert( is_numerically_one( detQ ) or is_numerically_one( -detQ ) );
                Assert( is_numerically_close( det, detQ * detR ) );
            }
            
            const auto Acof = CofactorMatrix( A );
            const auto Ainv = Inverse( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            LOG << Ainv << nl;
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;

            Assert( ( Ainv * A ).is_numerically_identity(), Ainv * A );
            Assert( ( A * Ainv ).is_numerically_identity(), A * Ainv );

            LOG << "with Newton-Schulz \n";
            auto Ainvp = Ainv; NewtonSchulz( A, Ainvp, 20 );

            LOG << A * Ainvp << nl;
            LOG << Ainvp * A << nl;

            auto R1 = ( A * Ainvp - IdentityMatrix(t) ).norm() / ( A * Ainv - IdentityMatrix(t) ).norm();
            auto R2 = ( Ainvp * A - IdentityMatrix(t) ).norm() / ( Ainv * A - IdentityMatrix(t) ).norm(); 
            
            LOG << "R1 = " << R1 << space << "R2 = " << R2 << nl;
            if( t >= 5 ) assert( R1 < 1.0 and R2 < 1.0 );

            assert( ( Ainvp * A ).is_numerically_identity() );
            assert( ( A * Ainvp ).is_numerically_identity() );
        }
    }
    

    
    {
        LOG << "7. Compare Determinants of random ill-conditioned matrices" << nl;
        for( int t = 2; t < 5; t++ )
        {
            
            DenseMatrix Q(2);
            Q.random_orthogonal_matrix();

            DenseMatrix M(2,2);
            M(0,0) = power_numerical( 10, t );
            M(0,1) = 0.;
            M(1,0) = 0.;
            M(1,1) = 1. / power_numerical( 10, t );

            DenseMatrix A = Transpose(Q) * M * Q;
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << t << " Determinant (default): " << det << nl;
            LOG << t << " Determinant (laplace): " << det_l << nl;
            LOG << t << " Determinant (gauss):   " << det_g << nl;
            LOG << t << " Determinant (bareiss): " << det_b << nl;

            assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );
            
            {
                const int dim = 2;
                DenseMatrix Q(dim,dim), R(dim,dim);
                QRFactorization( A, Q, R );
                Float detR = UpperTriangularDeterminant(R); Float detQ = Determinant(Q);
                LOG << "Determinant (QR): " << detR << nl;
                Assert( is_numerically_one( detQ ) or is_numerically_one( -detQ ) );
                Assert( is_numerically_close( det, detQ * detR ) );
            }
            
            const auto Acof = CofactorMatrix( A );
            const auto Ainv = Inverse( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;

            assert( ( Ainv * A ).is_numerically_identity() );
            assert( ( A * Ainv ).is_numerically_identity() );

            // LOG << "with Newton-Schulz \n";
            // auto Ainvp = Ainv; NewtonSchulz( A, Ainvp, 10 );

            // LOG << A * Ainvp << nl;
            // LOG << Ainvp * A << nl;

            // auto R1 = ( A * Ainvp - IdentityMatrix(t) ).norm() / ( A * Ainv - IdentityMatrix(t) ).norm();
            // auto R2 = ( Ainvp * A - IdentityMatrix(t) ).norm() / ( Ainv * A - IdentityMatrix(t) ).norm(); 
            
            // LOG << "R1 = " << R1 << space << "R2 = " << R2 << nl;
            // if( t >= 5 ) assert( R1 < 1.0 and R2 < 1.0 );

            // assert( ( Ainvp * A ).is_numerically_identity() );
            // assert( ( A * Ainvp ).is_numerically_identity() );
        }
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
    
    
}
