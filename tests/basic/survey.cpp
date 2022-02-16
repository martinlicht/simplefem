
#include <cstdio>
#include <cfloat>
#include <iostream>
#include "../../basic.hpp"

using namespace std;

int main()
{
    cout << "Survey of machine data (floating-point)" << endl;

    assert( machine_epsilon < 1e-10 );

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

    printf("\nProject-defined floating-point type\n");
    printf("Float decimal digits:   %ld\n",  (long int) std::numeric_limits<Float>::digits10         );
    printf("Float digits:           %ld\n",  (long int) std::numeric_limits<Float>::digits           );
    printf("Float min exponent:     %ld\n",  (long int) std::numeric_limits<Float>::min_exponent     );
    printf("Float max exponent:     %ld\n",  (long int) std::numeric_limits<Float>::max_exponent     );
    printf("Float min exponent 10:  %ld\n",  (long int) std::numeric_limits<Float>::min_exponent10   );
    printf("Float max exponent 10:  %ld\n",  (long int) std::numeric_limits<Float>::max_exponent10   );

    printf("Float denormalized min: %Le\n",  (long double) std::numeric_limits<Float>::denorm_min()  );
    printf("Float minimum:          %Le\n",  (long double) std::numeric_limits<Float>::min()         );
    printf("Float maximum:          %Le\n",  (long double) std::numeric_limits<Float>::max()         );
    printf("Float machine epsilon:  %Le\n",  (long double) std::numeric_limits<Float>::epsilon()     );

    printf("Float rounding error:   %Lf\n",  (long double) std::numeric_limits<Float>::round_error() );
    printf("Float has quiet NaN:     %d\n",  (int) std::numeric_limits<Float>::has_quiet_NaN         );
    printf("Float has signaling NaN: %d\n",  (int) std::numeric_limits<Float>::has_signaling_NaN     );
    // TODO: Show the same properties in the same order as above 

    return 0;
}
