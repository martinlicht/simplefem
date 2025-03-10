


#include <cstdio>

#include "../../basic.hpp"

// #include <cassert>
// #include <cstdlib>
// #include <cctype>
// #include <errno.h>
// #include <climits>
// #include <limits>

int main( int argc, char *argv[] ) 
{

    LOG << "Unit Test: string to integer conversion" << nl;
    
    if( argc == 2 ) {
    
        const char* endptr; 
        bool has_overflown;
        int value = string_to_integer( argv[1], &endptr, 10, has_overflown );

        bool invalid_end = ( *endptr != '\0' );
        
        if( has_overflown or invalid_end ) {

            if( has_overflown ) {
                fprintf( stderr, "Overflow detected: %s\n", argv[1] );
            }
        
            if( *endptr != '\0' ) {
                fprintf( stderr, "Invalid integer: %s\n", argv[1] );
            }

            return 1;
        }

        printf( "The integer value is: %d\n", value );
    
    } else {

        fprintf( stderr, "Usage: %s <integer>\n", argv[0]);
        return 0; // return 1;
        
    }

    

    {
        bool has_overflown;
        const char *endptr;
        int value;
    
        // Test 1: Simple positive number
        {
            const char* str = "123";
            value = string_to_integer( str, &endptr, 10, has_overflown);
            assert(value == 123);
            assert(*endptr == '\0');
            assert(!has_overflown);
        }   

        // Test 2: Simple negative number
        {
            const char* str = "-123";
            value = string_to_integer( str, &endptr, 10, has_overflown);
            assert(value == -123);
            assert(*endptr == '\0');
            assert(!has_overflown);
        }   

        // Test 3: Simple negative number with leading whitespace
        {
            const char* str = "  -456";
            value = string_to_integer(str, &endptr, 10, has_overflown);
            assert(value == -456);
            assert(*endptr == '\0');
            assert(!has_overflown);
        }   

        // Test 4: Number with trailing non-digit characters
        {
            const char* str = "789abc";
            value = string_to_integer( str, &endptr, 10, has_overflown);
            assert(value == 789);
            assert( endptr == str + 3 ); // The conversion stops at the first non-digit, so endptr should point to "abc"
            assert(!has_overflown);
        }   
        
        // Test 5: Hexadecimal conversion using base 16
        {
            const char* str = "1A";
            value = string_to_integer( str, &endptr, 16, has_overflown);
            assert(value == 26);
            assert(*endptr == '\0');
            assert(!has_overflown);
        }   
        
        // Test 6: Invalid digit for the given base
        {
            const char* str = "1G";
            value = string_to_integer( str, &endptr, 16, has_overflown);
            assert(value == 1);
            assert(*endptr == 'G'); // Conversion stops when encountering 'G'
            assert(!has_overflown);
        }   
        
        // Test 7: Positive overflow (assuming 32-bit int, INT_MAX = 2147483647)
        {
            const char* str = "2147483648";
            value = string_to_integer( str, &endptr, 10, has_overflown);
            assert(has_overflown);
            assert(*endptr == '\0');
            assert(value == std::numeric_limits<int>::max());
        }   
        
        // Test 8: Negative overflow (assuming INT_MIN = -2147483648)
        {
            const char* str = "-2147483649";
            value = string_to_integer( str, &endptr, 10, has_overflown);
            assert(has_overflown);
            assert(*endptr == '\0');
            assert(value == std::numeric_limits<int>::min());
        }   
        
        for( int base = 0; base <= 36; base++ ) {

            if( base == 1 ) continue;

            // Test 9a: Input with only whitespace should return 0
            {
                const char* str = "    ";
                value = string_to_integer( str, &endptr, 10, has_overflown);
                assert(value == 0);
                assert(*endptr == '\0');
                assert(!has_overflown);
            }   
            
            // Test 9b: Input with only a plus sign should return 0
            {
                const char* str = "+";
                value = string_to_integer( str , &endptr, 10, has_overflown);
                assert(value == 0);
                assert(*endptr == '\0'); // Since there's no digit, endptr should point to the string's end
                assert(!has_overflown);
            }   
            
            // Test 9c: Input with only a minus sign should return 0
            {
                const char* str = "-";
                value = string_to_integer( str , &endptr, 10, has_overflown);
                assert(value == 0);
                assert(*endptr == '\0'); // Since there's no digit, endptr should point to the string's end
                assert(!has_overflown);
            }   
            
            // Test 9d: Input with only empty string should return 0
            {
                const char* str = "";
                value = string_to_integer( str , &endptr, 10, has_overflown);
                assert(value == 0);
                assert(*endptr == '\0'); // Since there's no digit, endptr should point to the string's end
                assert(!has_overflown);
            }   
        
        }

    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    return 0;
}
