

#include "../../basic.hpp"
#include "../../dense/factorization.hpp"
#include "../../dense/functions.hpp"
// #include "../../dense/scalarfunctions.hpp"


// Function to compute the determinant of the Hilbert matrix
inline double hilbert_determinant(int n) {
    
    long double detinv = 1.;

    for (int k = 1; k <= n-1; k++) 
    {
        long double aux1 = factorial_numerical(2*k);
        long double aux2 = factorial_numerical(k)*factorial_numerical(k);
        long double binom = aux1/aux2;
        detinv *= (2*k+1) * binom*binom;
    }
    
    // LOG << "H1: " << 1./(double)detinv << std::endl;
    // The determinant is c_n / denominator
    return (double)1./detinv;
}




// Main function to test Vandermonde determinant computations
int main( int argc, char *argv[] )
{

    LOG << "Unit Test for Determinant computation" << nl;

    {
        LOG << nl << nl << "Vandermonde Matrices" << nl;

        // List of coefficient vectors for Vandermonde matrices
        std::vector<std::vector<Float>> coefficientVectors = {
            {7},
            {1, 2},
            {1, 2, 3},
            {1, 4, 9},
            {2, 3, 5, 7},
            {0.1, 0.2, 0.3, 0.4}
        };

        // Precomputed reference determinants
        std::vector<Float> referenceDeterminants = {
            1.0,          // empty product
            1.0,          // (2 - 1)
            2.0,          // (2 - 1)(3 - 1)(3 - 2)
            120.0,        // (4 - 1)(9 - 1)(9 - 4) = 3 * 8 * 5
            240.0,        // (3 - 2)(5 - 2)(5 - 3)(7 - 2)(7 - 3)(7 - 5) = 1 * 3 * 2 * 5 * 4 * 2 = 16*15 = 160+80 = 240
            12./1000000   // (2-1)(3-1)(4-1)(3-2)(4-2)(4-3) = 1 2 3 1 2 1 = 2 3 2 = 12 
        };

        assert( referenceDeterminants.size() == coefficientVectors.size() );

        // Test each Vandermonde matrix
        for( int i = 0; i < coefficientVectors.size(); ++i ) 
        {
            const auto& coeffs = coefficientVectors[i];
            Float referenceDet = referenceDeterminants[i];

            // Create and fill Vandermonde matrix
            DenseMatrix matrix( coeffs.size(), coeffs.size(), 1.);
            for( int r = 0; r < coeffs.size(); r++ )
            for( int c = 1; c < coeffs.size(); c++ )
            {
                matrix(r,c) = matrix(r,c-1) * coeffs[r];
            }

            // LOG << matrix << nl;

            // Compute determinant
            const int num_algorithms = 6;
            
            Float computedDet[num_algorithms];
            computedDet[0] = Determinant(matrix);
            computedDet[1] = Determinant_bareiss(matrix);
            computedDet[2] = Determinant_gauss(matrix);
            computedDet[3] = Determinant_laplaceexpansion(matrix);
            computedDet[4] = Determinant_ModifiedGramSchmidt(matrix);
            {
                DenseMatrix Q(matrix); DenseMatrix R(matrix);
                computedDet[5] = QRFactorization_via_Householder(matrix,Q,R);
                // LOG << matrix << nl << Q << nl << R << nl;
                // // LOG << matrix*Transpose(matrix) << nl << Q*Transpose(Q) << nl;
                // LOG << Transpose(matrix)*matrix << nl << Transpose(Q)*Q << nl;
                // LOG << machine_epsilon <<nl;
            }

            const char* names[num_algorithms];
            names[0] = "Determinant";
            names[1] = "Bareiss";
            names[2] = "Gauss";
            names[3] = "Laplace";
            names[4] = "Gram-Schmidt";
            names[5] = "Householder";

            // Print results
            for( int j = 0; j < num_algorithms; j++ )
            {
                LOGPRINTF("Test %i %12s: computed=%.17le reference=%.17le ratio=%.17le\n", 
                    i+1, names[j], referenceDet, computedDet[j], (computedDet[j]/referenceDet) );            
                assert( is_numerically_close( referenceDet, computedDet[j] ) );
            }
            LOG << nl;
        }    
    }
    
    
    {
        LOG << nl << nl << "Hilbert Matrices" << nl;

        for( int h = 0; h <= 8; h++ )
        {
            DenseMatrix matrix(h,h);
            for( int i = 0; i < h; i++ ) 
            for( int j = 0; j < h; j++ ) 
                matrix(i,j) = 1. / ( i + j + 1 );

            Float referenceDet = hilbert_determinant(h);

            // LOG << matrix << nl;

            // Compute determinant
            const int num_algorithms = 6;
            
            Float computedDet[num_algorithms];
            computedDet[0] = Determinant(matrix);
            computedDet[1] = Determinant_bareiss(matrix);
            computedDet[2] = Determinant_gauss(matrix);
            computedDet[3] = Determinant_laplaceexpansion(matrix);
            computedDet[4] = Determinant_ModifiedGramSchmidt(matrix);
            {
                DenseMatrix Q(matrix); DenseMatrix R(matrix);
                computedDet[5] = QRFactorization_via_Householder(matrix,Q,R);
                LOG << matrix << nl << Q << nl << R << nl;
                // // LOG << matrix*Transpose(matrix) << nl << Q*Transpose(Q) << nl;
                // LOG << Transpose(matrix)*matrix << nl << Transpose(Q)*Q << nl;
                // LOG << machine_epsilon <<nl;
            }
            
            const char* names[num_algorithms];
            names[0] = "Determinant";
            names[1] = "Bareiss";
            names[2] = "Gauss";
            names[3] = "Laplace";
            names[4] = "Gram-Schmidt";
            names[5] = "Householder";

            // Print results
            for( int j = 0; j < num_algorithms; j++ )
            {
                LOGPRINTF("Test %i %12s: reference=%.17le computed=%.17le ratio=%.17le\n", 
                    h, names[j], referenceDet, computedDet[j], (computedDet[j]/(double)referenceDet) );            
                assert( is_numerically_close( referenceDet, computedDet[j] ) );
            }
            LOG << nl;

        }

    }

    

    {
        LOG << nl << nl << "Random orthogonal matrices" << nl;

        for( int h = 0; h <= 8; h++ )
        {
            DenseMatrix matrix(h,h);
            matrix.random_orthogonal_matrix();
            Float referenceDet = 1.;

            matrix.randommatrix();
            for( int i = 0; i < 1; i++ )
            {
                DenseMatrix Q(matrix); DenseMatrix R(matrix);
                QRFactorization_via_Householder(matrix,Q,R);
                matrix = Q;
            }
            
            // LOG << matrix << nl;

            // Compute determinant
            const int num_algorithms = 6;
            
            Float computedDet[num_algorithms];
            computedDet[0] = Determinant(matrix);
            computedDet[1] = Determinant_bareiss(matrix);
            computedDet[2] = Determinant_gauss(matrix);
            computedDet[3] = Determinant_laplaceexpansion(matrix);
            computedDet[4] = Determinant_ModifiedGramSchmidt(matrix);
            {
                DenseMatrix Q(matrix); DenseMatrix R(matrix);
                computedDet[5] = QRFactorization_via_Householder(matrix,Q,R);
                // LOG << matrix << nl << Q << nl << R << nl;
                // // LOG << matrix*Transpose(matrix) << nl << Q*Transpose(Q) << nl;
                // LOG << Transpose(matrix)*matrix << nl << Transpose(Q)*Q << nl;
                // LOG << machine_epsilon <<nl;
            }
            
            
            const char* names[num_algorithms];
            names[0] = "Determinant";
            names[1] = "Bareiss";
            names[2] = "Gauss";
            names[3] = "Laplace";
            names[4] = "Gram-Schmidt";
            names[5] = "Householder";

            // Print results
            for( int j = 0; j < num_algorithms; j++ )
            {
                computedDet[j] = absolute( computedDet[j] );
                LOGPRINTF("Test %i %12s: reference=%.17le computed=%.17le ratio=%.17le\n", 
                    h, names[j], referenceDet, computedDet[j], (computedDet[j]/referenceDet) );            
                assert( is_numerically_close( referenceDet, computedDet[j] ) );
            }
            LOG << nl;

        }

    }


    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    return 0;
}
