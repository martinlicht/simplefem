#include <iostream>
#include <cassert>

#include "../../basic.hpp"
#include "../../utility/random.hpp"

// A small tolerance for floating-point comparisons

int main( int argc, char *argv[] )
{

    LOG << "Unit Test for randomness" << nl;

    {
        LOG << "testing random integer generation...\n";
        seed_random_integer();
        
        unsigned int max_value = random_integer_maximum();
        
        // Ensure `random_integer_maximum()` returns a non-zero maximum
        
        // Ensure `random_integer` returns values within a reasonable range
        assert( max_value > 0 && "random_integer_maximum() should be greater than 0" );

        for( int t = 0; t < 10; t++ ) {
            unsigned int value = random_integer();
            unsigned int max   = random_integer_maximum();
            assert( max == max_value );
            assert( value <= max_value && "random_integer() out of range" );
        }
        
    }

    {
        std::cout << "Test coin flip...\n";
        seed_random_integer();
        
        // Test flipping a coin with a 50% probability
        unsigned int result = flip_coin();
        assert( (result == 0 || result == 1) && "flip_coin() must return 0 or 1" );

        // Test flipping a coin with a 100% probability of zero
        result = flip_coin(1.0);
        assert( result == 0 && "flip_coin(1.0) should always return 0" );

        // Test flipping a coin with a 0% probability of zero
        result = flip_coin(0.0);
        assert( result == 1 && "flip_coin(0.0) should always return 1" );

        // test whether coin is fair
        int count[2] = { 0, 0 };
        const int M = 1 << 10;
        for( int t = 0; t < M; t++ )
        {
            const int c = flip_coin();
            assert( c == 0 or c == 1 );
            count[c]++;
        }

        const Float TOLERANCE = 1e-1;

        double share0 = count[0]/(double)M; 
        double share1 = count[1]/(double)M; 
        LOG << share0 << space << share1 << nl;
        assert( std::abs(share0-0.5) < TOLERANCE );
        assert( std::abs(share1-0.5) < TOLERANCE );
    }



    {
        std::cout << "test random uniform...\n";

        seed_random_integer();
        
        // Generate a random uniform number
        for( int t = 0; t < 1000; t++ )
        {
            Float value = random_uniform();
            assert( value >= 0.0 && value <= 1.0 && "random_uniform() must return values in [0, 1]" );
        }
        
    }

    {
        std::cout << "test_gaussrand...\n";
        seed_random_integer();
        
        // Generate a random Gaussian number
        for( int t = 0; t < 1000; t++ )
        {
            Float value = gaussrand();

            const Float TOLERANCE = 1e+6;

            // Basic sanity check (this doesn't verify the distribution, but it's a start)
            assert( value == value && "gaussrand() should not return NaN" ); // Ensure no NaN
            assert( std::abs(value) < TOLERANCE && "gaussrand() value too extreme" );
        }
        
    }
        
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
    return 0;
}
