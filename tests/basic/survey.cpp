
#include <cstdio>
#include <cfloat>
#include <iostream>
#include "../../basic.hpp"

using namespace std;

int main()
{
    cout << "Survey of machine data" << endl;

    /* output floating-point parameters */

    printf("    \nFloating-point mode flags\n");
    printf("    FLT_ROUNDS      : %ld\n", (long) FLT_ROUNDS      );
    printf("    FLT_EVAL_METHOD : %ld\n", (long) FLT_EVAL_METHOD );
    printf("    FLT_RADIX       : %ld\n", (long) FLT_RADIX       );
    printf("    DECIMAL_DIG     : %ld\n", (long) DECIMAL_DIG     );

    printf("    \nFloating-point properties\n");
    printf("    FLT/DBL/LDBL_DECIMAL_DIG : %ld\t%ld\t%ld\n", (long) FLT_DECIMAL_DIG, (long) DBL_DECIMAL_DIG, (long) LDBL_DECIMAL_DIG );
    printf("    FLT/DBL/LDBL_DIG         : %ld\t%ld\t%ld\n", (long) FLT_DIG,         (long) DBL_DIG,         (long) LDBL_DIG         );
    printf("    FLT/DBL/LDBL_MANT_DIG    : %ld\t%ld\t%ld\n", (long) FLT_MANT_DIG,    (long) DBL_MANT_DIG,    (long) LDBL_MANT_DIG    );
    printf("    FLT/DBL/LDBL_MIN_EXP     : %ld\t%ld\t%ld\n", (long) FLT_MIN_EXP,     (long) DBL_MIN_EXP,     (long) LDBL_MIN_EXP     );
    printf("    FLT/DBL/LDBL_MAX_EXP     : %ld\t%ld\t%ld\n", (long) FLT_MAX_EXP,     (long) DBL_MAX_EXP,     (long) LDBL_MAX_EXP     );
    printf("    FLT/DBL/LDBL_MIN_10_EXP  : %ld\t%ld\t%ld\n", (long) FLT_MIN_10_EXP,  (long) DBL_MIN_10_EXP,  (long) LDBL_MIN_10_EXP  );
    printf("    FLT/DBL/LDBL_MAX_10_EXP  : %ld\t%ld\t%ld\n", (long) FLT_MAX_10_EXP,  (long) DBL_MAX_10_EXP,  (long) LDBL_MAX_10_EXP  );
    printf("    FLT/DBL/LDBL_HAS_SUBNORM : %ld\t%ld\t%ld\n", (long) FLT_HAS_SUBNORM, (long) DBL_HAS_SUBNORM, (long) LDBL_HAS_SUBNORM );

    printf("    \nFloating-point Minima and Maxima\n");
    printf("    FLT/DBL/LDBL_TRUE_MIN : %Le\t%Le\t%Le\n", (long double) FLT_TRUE_MIN, (long double) DBL_TRUE_MIN, (long double) LDBL_TRUE_MIN );
    printf("    FLT/DBL/LDBL_MIN      : %Le\t%Le\t%Le\n", (long double) FLT_MIN,      (long double) DBL_MIN,      (long double) LDBL_MIN      );
    printf("    FLT/DBL/LDBL_MAX      : %Le\t%Le\t%Le\n", (long double) FLT_MAX,      (long double) DBL_MAX,      (long double) LDBL_MAX      );
    printf("    FLT/DBL/LDBL_EPSILON  : %Le\t%Le\t%Le\n", (long double) FLT_EPSILON,  (long double) DBL_EPSILON,  (long double) LDBL_EPSILON  );

    printf("\nProject-defined floating-point type");
    printf("Float max exponent:     %ld\n",  (long int) std::numeric_limits<Float>::max_exponent     );
    printf("Float min exponent:     %ld\n",  (long int) std::numeric_limits<Float>::min_exponent     );
    printf("Float max exponent 10:  %ld\n",  (long int) std::numeric_limits<Float>::max_exponent10   );
    printf("Float min exponent 10:  %ld\n",  (long int) std::numeric_limits<Float>::min_exponent10   );
    printf("Float denormalized min: %Le\n",  (long double) std::numeric_limits<Float>::denorm_min()  );
    printf("Float rounding error:   %Lf\n",  (long double) std::numeric_limits<Float>::round_error() );
    printf("Float machine epsilon:  %Le\n",  (long double) std::numeric_limits<Float>::epsilon()     );
    printf("Float minimum:          %Le\n",  (long double) std::numeric_limits<Float>::min()         );
    printf("Float maximum:          %Le\n",  (long double) std::numeric_limits<Float>::max()         );
    printf("Float has quiet NaN:     %d\n",  (int) std::numeric_limits<Float>::has_quiet_NaN         );
    printf("Float has signaling NaN: %d\n",  (int) std::numeric_limits<Float>::has_signaling_NaN     );
    // TODO: Show the same properties in the same order as above 

    /* survey factorials */

    std::cout << "\nLargest number whose factorial fits into data type (signed)" << std::endl;
    std::cout << "    case signed char       : " << largest_factorial_base<       signed char>() << std::endl;
    std::cout << "    case signed short      : " << largest_factorial_base<      signed short>() << std::endl;
    std::cout << "    case signed int        : " << largest_factorial_base<        signed int>() << std::endl;
    std::cout << "    case signed long       : " << largest_factorial_base<       signed long>() << std::endl;
    std::cout << "    case signed long long  : " << largest_factorial_base<  signed long long>() << std::endl;

    std::cout << "\nLargest number whose factorial fits into data type (unsigned)" << std::endl;
    std::cout << "    case unsigned char     : " << largest_factorial_base<     unsigned char>() << std::endl;
    std::cout << "    case unsigned short    : " << largest_factorial_base<    unsigned short>() << std::endl;
    std::cout << "    case unsigned int      : " << largest_factorial_base<      unsigned int>() << std::endl;
    std::cout << "    case unsigned long     : " << largest_factorial_base<     unsigned long>() << std::endl;
    std::cout << "    case unsigned long long: " << largest_factorial_base<unsigned long long>() << std::endl;

    for( int n = 0; n < 50; n++ ) {
        auto x = largest_factorial_base_AUX( n, 2 );
        //printf( "Largest k such that k! <= %d : %lld\n", n, (long long int)x );
        if(  0 <= x && x <=   1 ) Assert( x == 1 );
        if(  2 <= x && x <=   4 ) Assert( x == 2 );
        if(  5 <= x && x <=  23 ) Assert( x == 3 );
        if( 24 <= x && x <= 124 ) Assert( x == 4 );
    }

    /* TODO: check factorials, binomials, */

    return 0;
}
