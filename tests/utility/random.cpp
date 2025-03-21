
#include <cmath>

#include "../../base/include.hpp"
#include "../../utility/random.hpp"

// A small tolerance for floating-point comparisons

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: randomness" << nl;

    {
        LOG << "testing random integer generation...\n";
        seed_random_integer();
        
        unsigned int max_value =  get_random_integer_modulo();
        
        // Ensure ` get_random_integer_modulo()` returns a non-zero maximum
        assert( max_value > 0 );

        // Ensure `random_integer` returns values within a reasonable range        
        for( int t = 0; t < 10; t++ ) {
            unsigned int value = random_integer();
            unsigned int max   =  get_random_integer_modulo();
            assert( max == max_value );
            assert( value <= max_value );
        }
    }

    {
        LOG << "Test coin flip...\n";
        seed_random_integer();
        
        unsigned int result = flip_coin();

        // Test flipping a coin with a 50% probability
        
        for( int i = 0; i < 100; i++ ) {
            
            // The result must be either 0 or 1
            result = flip_coin();
            assert( result == 0 || result == 1 );

            // Test flipping a coin with a 100% probability of zero
            result = flip_coin(1.0);
            assert( result == 0 );

            // Test flipping a coin with a 0% probability of zero
            result = flip_coin(0.0);
            assert( result == 1 );

        }

        // test whether the coin is fair
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
        LOG << "test random uniform...\n";

        seed_random_integer();
        
        // Generate a random uniform number
        for( int t = 0; t < 1000; t++ )
        {
            Float value = random_uniform();
            assert( value >= 0.0 && value <= 1.0 && "random_uniform() must return values in [0, 1]" );
        }
    }

    {
        LOG << "test_gaussian_variable...\n";
        seed_random_integer();
        
        // Generate a random Gaussian number
        for( int t = 0; t < 1000; t++ )
        {
            Float value = gaussian_variable();

            const Float TOLERANCE = 1e+6;

            // Basic sanity check (not verifying the distribution but it is a start)
            assert( value == value && "gaussian_variable() should not return NaN" ); // Ensure no NaN
            assert( std::abs(value) < TOLERANCE && "gaussian_variable() value too extreme" );
        }
    }
        
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
    return 0;
}
