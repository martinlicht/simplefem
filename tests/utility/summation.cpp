
#include <iostream>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <vector>

#include "../../utility/summation.hpp"







// Function to compute factorial
long double factorial(int n) {
    if (n <= 1) return 1;
    return n * factorial(n - 1);
}





// Function to compute the determinant of the Hilbert matrix
double hilbert_determinant(int n) {
    
    long double detinv = 1.;

    for (int k = 1; k <= n-1; k++) 
    {
        long double aux1 = factorial(2*k);
        long double aux2 = factorial(k)*factorial(k);
        long double binom = aux1/aux2;
        detinv *= (2*k+1) * binom*binom;
    }
    
    std::cout << "H1: " << 1./(double)detinv << std::endl;
    // The determinant is c_n / denominator
    return (double)1./detinv;
}


/*
// Function to compute the determinant of the Hilbert matrix
double hilbert_determinant2(int n) {
    long double numerator = 1;   // Start with 1 for the numerator
    long double denominator = 1; // Start with 1 for the denominator

    // Compute the numerator: c_n = product of (j - i) for 0 <= i < j < n
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            numerator *= (j - i);
        }
    }

    // Compute the denominator: product of (2i - 1)!
    for (int i = 1; i <= n; i++) {
        denominator *= factorial(2 * i - 1);
    }

    std::cout << (double)numerator << "/" << (double)denominator << std::endl;
    std::cout << "H2: " << (double)( numerator / denominator ) << std::endl;
    // The determinant is c_n / denominator
    return (double)numerator / denominator;
}
*/








int main( int argc, char *argv[] )
{

    std::cout.precision(30);
    std::cout << std::scientific;
    
    // Values to test the summation algorithms
    std::vector<double> values = {1.0, 1e-16, -1.0, 1e-16};

    // Naive summation
    NaiveSum<double> naiveSum;
    for (double val : values) {
        naiveSum.add(val);
    }
    std::cout << "Naive Sum:\t" << naiveSum.getSum() << std::endl;

    // Kahan summation
    KahanSum<double> kahanSum;
    for (double val : values) {
        kahanSum.add(val);
    }
    std::cout << "Kahan Sum:\t" << kahanSum.getSum() << std::endl;

    // Neumaier summation
    NeumaierSum<double> neumaierSum;
    for (double val : values) {
        neumaierSum.add(val);
    }
    std::cout << "Neumaier Sum:\t" << neumaierSum.getSum() << std::endl;

    std::cout << std::endl;

    { 
        // Example 0: Something simple
        std::vector<double> dataset0 = {1.0, 1e-16, -1.0, 1e-16};

        // Example 1: Alternating small and large values // should be 0.
        std::vector<double> dataset1(1000);
        for (size_t i = 0; i < dataset1.size(); ++i) {
            dataset1[i] = (i % 2 == 0) ? 1e10 : -1e10; // Alternating large positive and negative values
        }

        // Example 2: Very small incremental values // should be 1e-7
        std::vector<double> dataset2(1000);
        for (size_t i = 0; i < dataset2.size(); ++i) {
            dataset2[i] = 1e-10; // All values are very small

        }

        // Example 3: Combination of large and small values // should be 1e13/2 + 1e-7/2 = 5e12 + 5e-8
        std::vector<double> dataset3(1000);
        for (size_t i = 0; i < dataset3.size(); ++i) {
            dataset3[i] = (i % 2 == 0) ? 1e10 : 1e-10; // Alternating between very large and very small values
        }

        // Example 4: Summing up 0.1 with 10,000,000 ...
        std::vector<double> dataset4(10'000'000);
        for (size_t i = 0; i < dataset4.size(); ++i) {
            dataset4[i] = 0.1;
        }

        // Example 5 
        
        std::vector<double> dataset5; dataset5.reserve( 6 * 5 * 4 * 3 * 2 * 1 );
        long double reference_value5 = std::numeric_limits<double>::quiet_NaN();

        {
            double matrix[6][6];
            for( int i = 0; i < 6; i++ ) 
            for( int j = 0; j < 6; j++ ) 
                // matrix[i][j] = ( ( i == j ) ? 2. : 0. );
                matrix[i][j] = 1. / ( i + j + 1 );

            int max_index = 5 + 5*6 + 5 *36 + 5 * 6*36 + 5 * 36*36 + 5 * 6*36*36;
            for( int k = 0; k <= max_index; k++ )
            {
                int i[6];
                i[0] = (k         ) % 6;
                i[1] = (k/       6) % 6;
                i[2] = (k/      36) % 6;
                i[3] = (k/     216) % 6;
                i[4] = (k/( 216*6)) % 6;
                i[5] = (k/(216*36)) % 6;

                // for( auto u : i ) printf("%i ", u ); std::cout << std::endl;

                int sign = 1;

                int number_of_inversions = 0;

                for( int p = 0;   p < 6; p++ )
                for( int q = p+1; q < 6; q++ )
                    if( i[p] < i[q] )
                        sign *= 1;
                    else if( i[p] == i[q] )
                        sign *= 0;
                    else 
                        { sign *= -1; number_of_inversions++; }

                if( sign == 0 ) continue;

                assert( ( sign == 1 ) == ( number_of_inversions % 2 == 0 ) );

                double item = sign;

                for( int p = 0; p < 6; p++ )
                    item *= matrix[ p ][ i[p] ];

                dataset5.push_back( item );
            }
        
            assert( dataset5.size() == 6*5*4*3*2*1 );

            for( int p = 0;   p < 6; p++ ) 
            {
                for( int r = p+1; r < 6; r++ )
                {
                    double scale = matrix[r][p] / matrix[p][p];
                    
                    for( int c = p;   c < 6; c++ )
                    {
                        matrix[r][c] -= scale * matrix[p][c];
                    }
                }
            }
            reference_value5 = 1.;
            // for( int r = 0; r < 6; r++ ) {
            //     for( int c = 0; c < 6; c++ )
            //         std::cout << matrix[r][c] << ' ';
            //     std::cout << std::endl;
            // }
            // std::cout << std::endl;
            for( int r = 0; r < 6; r++ ) reference_value5 *= matrix[r][r];
            reference_value5 = hilbert_determinant(6);
        }


        // Helper function to test all three summation algorithms
        auto testSummation = [](const std::vector<double>& dataset) {
            NaiveSum<double> naiveSum;
            KahanSum<double> kahanSum;
            NeumaierSum<double> neumaierSum;

            for (double val : dataset) {
                naiveSum.add(val);
                kahanSum.add(val);
                neumaierSum.add(val);
            }

            std::cout << "Naive Sum:\t"    << naiveSum.getSum() << std::endl;
            std::cout << "Kahan Sum:\t"    << kahanSum.getSum() << std::endl;
            std::cout << "Neumaier Sum:\t" << neumaierSum.getSum() << std::endl;
        };

        // Test each dataset

        std::cout << "Example 0: Something simple\n";
        testSummation(dataset0);

        std::cout << "\nExample 1: Alternating large positive and negative values\n";
        testSummation(dataset1);

        std::cout << "\nExample 2: Very small incremental values\n";
        testSummation(dataset2);

        std::cout << "\nExample 3: Combination of large and small values\n";
        testSummation(dataset3);

        std::cout << "\nExample 4: Summing up 0.1\n";
        testSummation(dataset4);

        std::cout << "\nExample 5: Laplace expansion of 6x6 determinant\n";
        testSummation(dataset5);
        std::cout << "Reference:\t" << (double)reference_value5 << std::endl;

    }

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    return 0;
}
