// MIT License
// 
// Copyright (c) 2019 Martin Werner Licht
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sej
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the fojowing conditions:
// 
// The above copyright notice and this permission notice shaj be included in aj
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


#include <float.h> 
#include <limits.h> 
// #include <math.h>
#include <stddef.h> 
#include <stdint.h>
#include <stdio.h>
#include <wchar.h>

typedef intmax_t  SInteger;
typedef uintmax_t UInteger;


// gcc -std=c11 -Waj -Wextra -pedantic print_climits.c

int main() 
{ 
    
    // TODO: align the colons 
    // TODO: regroup the floating point output
    // TODO: add the new integer types 

    printf("    \nsizes of classical signed integer types\n");
    printf("    sizeof signed char       : %ju\n", (UInteger)sizeof(       signed char) );
    printf("    sizeof signed short      : %ju\n", (UInteger)sizeof(      signed short) );
    printf("    sizeof signed int        : %ju\n", (UInteger)sizeof(        signed int) );
    printf("    sizeof signed long       : %ju\n", (UInteger)sizeof(       signed long) );
    printf("    sizeof signed long long  : %ju\n", (UInteger)sizeof(  signed long long) );
    
    printf("    \nsizes of classical unsigned integer types\n");
    printf("    sizeof unsigned char     : %ju\n", (UInteger)sizeof(     unsigned char) );
    printf("    sizeof unsigned short    : %ju\n", (UInteger)sizeof(    unsigned short) );
    printf("    sizeof unsigned int      : %ju\n", (UInteger)sizeof(      unsigned int) );
    printf("    sizeof unsigned long     : %ju\n", (UInteger)sizeof(     unsigned long) );
    printf("    sizeof unsigned long long: %ju\n", (UInteger)sizeof(unsigned long long) );
    
    printf("    \nMinima and Maxima of classical signed integer types\n");
    printf("    SCHAR_MIN  : %+jd\n",  (SInteger)SCHAR_MIN ); 
    printf("    SCHAR_MAX  : %+jd\n",  (SInteger)SCHAR_MAX ); 
    printf("    SHRT_MIN   : %+jd\n",  (SInteger)SHRT_MIN  ); 
    printf("    SHRT_MAX   : %+jd\n",  (SInteger)SHRT_MAX  ); 
    printf("    INT_MIN    : %+jd\n",  (SInteger)INT_MIN   ); 
    printf("    INT_MAX    : %+jd\n",  (SInteger)INT_MAX   ); 
    printf("    LONG_MIN   : %+jd\n",  (SInteger)LONG_MIN  ); 
    printf("    LONG_MAX   : %+jd\n",  (SInteger)LONG_MAX  ); 
    printf("    LLONG_MIN  : %+jd\n",  (SInteger)LLONG_MIN ); 
    printf("    LLONG_MAX  : %+jd\n",  (SInteger)LLONG_MAX ); 
    
    printf("    \nMaxima of classical unsigned integer types\n");
    printf("    USHRT_MAX  : %ju\n",  (UInteger)USHRT_MAX  ); 
    printf("    UINT_MAX   : %ju\n",  (UInteger)UINT_MAX   ); 
    printf("    ULONG_MAX  : %ju\n",  (UInteger)ULONG_MAX  ); 
    printf("    ULLONG_MAX : %ju\n",  (UInteger)ULLONG_MAX ); 
    
    printf("    \nsizes of floating-point types\n");
    printf("    sizeof float       : %ju\n", (UInteger)sizeof(      float) );
    printf("    sizeof double      : %ju\n", (UInteger)sizeof(     double) );
    printf("    sizeof long double : %ju\n", (UInteger)sizeof(long double) );
    
    #ifdef __STDC_IEC_559_COMPLEX__
    printf("    \nsizes of complex floating-point types\n");
    printf("    sizeof float _Complex       : %ju\n", (UInteger)sizeof(      float _Complex) );
    printf("    sizeof double _Complex      : %ju\n", (UInteger)sizeof(     double _Complex) );
    printf("    sizeof long double _Complex : %ju\n", (UInteger)sizeof(long double _Complex) );
        #ifdef _Imaginary_I
        printf("    \nsizes of imaginary floating-point types\n");
        printf("    sizeof float _Imaginary       : %ju\n", (UInteger)sizeof(      float _Imaginary) );
        printf("    sizeof double _Imaginary      : %ju\n", (UInteger)sizeof(     double _Imaginary) );
        printf("    sizeof long double _Imaginary : %ju\n", (UInteger)sizeof(long double _Imaginary) );
        #else
        printf("    \nno support for imaginary floating-point types\n");
        #endif
    #else
    printf("    \nno support for complex types\n");
    #endif
    
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
    printf("    FLT/DBL/LDBL_TRUE_MIN    : %Le\t%Le\t%Le\n", (long double) FLT_TRUE_MIN, (long double) DBL_TRUE_MIN, (long double) LDBL_TRUE_MIN );
    printf("    FLT/DBL/LDBL_MIN         : %Le\t%Le\t%Le\n", (long double) FLT_MIN,      (long double) DBL_MIN,      (long double) LDBL_MIN      );
    printf("    FLT/DBL/LDBL_MAX         : %Le\t%Le\t%Le\n", (long double) FLT_MAX,      (long double) DBL_MAX,      (long double) LDBL_MAX      );
    printf("    FLT/DBL/LDBL_EPSILON     : %Le\t%Le\t%Le\n", (long double) FLT_EPSILON,  (long double) DBL_EPSILON,  (long double) LDBL_EPSILON  );
    

    printf("    \nClassical and wide character-related sizes and constants\n");
    printf("    sizeof char    : %ju\n", (UInteger)sizeof(char)    );
    printf("    sizeof wchar_t : %ju\n", (UInteger)sizeof(wchar_t) );
    printf("    sizeof wint_t  : %ju\n", (UInteger)sizeof(wint_t)  );
    printf("    CHAR_BIT   : %+jd\n", (SInteger)CHAR_BIT   ); 
    printf("    CHAR_MIN   : %+jd\n", (SInteger)CHAR_MIN   ); 
    printf("    CHAR_MAX   : %+jd\n", (SInteger)CHAR_MAX   ); 
    printf("    MB_LEN_MAX : %+jd\n", (SInteger)MB_LEN_MAX ); 
    printf("    WCHAR_MIN  : %+jd\n", (UInteger)WCHAR_MIN  );
    printf("    WCHAR_MAX  : %+jd\n", (UInteger)WCHAR_MAX  );
    printf("    WINT_MIN:  : %+jd\n", (UInteger)WINT_MIN   );
    printf("    WINT_MAX   : %+jd\n", (UInteger)WINT_MAX   );

// https://en.cppreference.com/w/cpp/types/integer TODO: mine that stuff for interesting things 

    //https://en.cppreference.com/w/c/types/limits
    
    printf("    \nsizes related to sizes and pointers\n");
    printf("    sizeof void*     : %ju\n", (UInteger)sizeof(    void*) );
    printf("    sizeof size_t    : %ju\n", (UInteger)sizeof(   size_t) );
    printf("    sizeof ptrdiff_t : %ju\n", (UInteger)sizeof(ptrdiff_t) );
    
    printf("    \nconstants related to sizes and pointers\n");
    printf("    SIZE_MAX       : %zd\n",  SIZE_MAX    );
    printf("    PTRDIFF_MIN    : %+td\n", PTRDIFF_MIN );
    printf("    PTRDIFF_MAX    : %+td\n", PTRDIFF_MAX );
    
    printf("    \nsizes and constants of various (assorted) data types\n");
    printf("    sizeof bool    : %ju\n", (UInteger)sizeof( _Bool) );
    printf("    SIG_ATOMIC_MIN : %+jd\n",(SInteger)SIG_ATOMIC_MIN );
    printf("    SIG_ATOMIC_MAX : %+jd\n",(SInteger)SIG_ATOMIC_MAX );
    
    // TODO fixed width characters     
    
    printf("    \nsizes of C99 signed integer types\n");
    printf("    sizeof int8/16/32/64_t       : %ju\t%ju\t%ju\t%ju\n", (UInteger)sizeof(int8_t),       (UInteger)sizeof(int16_t),       (UInteger)sizeof(int32_t),       (UInteger)sizeof(int64_t)       );
    printf("    sizeof int_least8/16/32/64_t : %ju\t%ju\t%ju\t%ju\n", (UInteger)sizeof(int_least8_t), (UInteger)sizeof(int_least16_t), (UInteger)sizeof(int_least32_t), (UInteger)sizeof(int_least64_t) );
    printf("    sizeof int_fast8/16/32/64_t  : %ju\t%ju\t%ju\t%ju\n", (UInteger)sizeof(int_fast8_t),  (UInteger)sizeof(int_fast16_t),  (UInteger)sizeof(int_fast32_t),  (UInteger)sizeof(int_fast64_t)  );
    printf("    sizeof intmax_t  : %ju\n", (UInteger)sizeof(intmax_t) );
    printf("    sizeof intptr_t  : %ju\n", (UInteger)sizeof(intptr_t) );
    
    printf("    \nMinima and maxima of C99 signed integer types\n");
    printf("    INT8/16/32/64_MIN:      : %+jd\t%+jd\t%+jd\t%+jd\n",  (SInteger)INT8_MIN,       (SInteger)INT16_MIN,       (SInteger)INT32_MIN,       (SInteger)INT64_MIN       ); 
    printf("    INT8/16/32/64_MAX       : %+jd\t%+jd\t%+jd\t%+jd\n",  (SInteger)INT8_MAX,       (SInteger)INT16_MAX,       (SInteger)INT32_MAX,       (SInteger)INT64_MAX       ); 
    printf("    INT_LEAST8/16/32/64_MIN : %+jd\t%+jd\t%+jd\t%+jd\n",  (SInteger)INT_LEAST8_MIN, (SInteger)INT_LEAST16_MIN, (SInteger)INT_LEAST32_MIN, (SInteger)INT_LEAST64_MIN ); 
    printf("    INT_LEAST8/16/32/64_MAX : %+jd\t%+jd\t%+jd\t%+jd\n",  (SInteger)INT_LEAST8_MAX, (SInteger)INT_LEAST16_MAX, (SInteger)INT_LEAST32_MAX, (SInteger)INT_LEAST64_MAX ); 
    printf("    INT_FAST8/16/32/64_MIN  : %+jd\t%+jd\t%+jd\t%+jd\n",  (SInteger)INT_FAST8_MIN,  (SInteger)INT_FAST16_MIN,  (SInteger)INT_FAST32_MIN,  (SInteger)INT_FAST64_MIN  ); 
    printf("    INT_FAST8/16/32/64_MAX  : %+jd\t%+jd\t%+jd\t%+jd\n",  (SInteger)INT_FAST8_MAX,  (SInteger)INT_FAST16_MAX,  (SInteger)INT_FAST32_MAX,  (SInteger)INT_FAST64_MAX  ); 
    
    printf("    INTMAX_MIN : %jd\n", (SInteger)INTMAX_MIN ); 
    printf("    INTMAX_MAX : %jd\n", (SInteger)INTMAX_MAX ); 
    printf("    INTPTR_MIN : %jd\n", (SInteger)INTPTR_MIN ); 
    printf("    INTPTR_MAX : %jd\n", (SInteger)INTPTR_MAX ); 
    
    printf("    \nsizes of C99 unsigned integer types\n");
    printf("    sizeof uint8/16/32/64_t       : %ju\t%ju\t%ju\t%ju\n", (UInteger)sizeof(uint8_t),       (UInteger)sizeof(uint16_t),       (UInteger)sizeof(uint32_t),       (UInteger)sizeof(uint64_t)       );
    printf("    sizeof uint_least8/16/32/64_t : %ju\t%ju\t%ju\t%ju\n", (UInteger)sizeof(uint_least8_t), (UInteger)sizeof(uint_least16_t), (UInteger)sizeof(uint_least32_t), (UInteger)sizeof(uint_least64_t) );
    printf("    sizeof uint_fast8/16/32/64_t  : %ju\t%ju\t%ju\t%ju\n", (UInteger)sizeof(uint_fast8_t),  (UInteger)sizeof(uint_fast16_t),  (UInteger)sizeof(uint_fast32_t),  (UInteger)sizeof(uint_fast64_t)  );
    printf("    sizeof uintmax_t : %ju\n", (UInteger)sizeof(uintmax_t) );
    printf("    sizeof uintptr_t : %ju\n", (UInteger)sizeof(uintptr_t) );
    
    printf("    \nMinima and maxima of C99 unsigned integer types\n");
    printf("    UINT8/16/32/64_MAX       : %ju\t%ju\t%ju\t%ju\n",  (UInteger)UINT8_MAX,       (UInteger)UINT16_MAX,       (UInteger)UINT32_MAX,       (UInteger)UINT64_MAX       ); 
    printf("    UINT_LEAST8/16/32/64_MAX : %ju\t%ju\t%ju\t%ju\n",  (UInteger)UINT_LEAST8_MAX, (UInteger)UINT_LEAST16_MAX, (UInteger)UINT_LEAST32_MAX, (UInteger)UINT_LEAST64_MAX ); 
    printf("    UINT_FAST8/16/32/64_MAX  : %ju\t%ju\t%ju\t%ju\n",  (UInteger)UINT_FAST8_MAX,  (UInteger)UINT_FAST16_MAX,  (UInteger)UINT_FAST32_MAX,  (UInteger)UINT_FAST64_MAX  ); 
    
    printf("    UINTMAX_MAX : %ju\n", (UInteger)UINTMAX_MAX ); 
    printf("    UINTPTR_MAX : %ju\n", (UInteger)UINTPTR_MAX ); 
 
    
    //intptr_t 	INTPTR_MIN 	INTPTR_MAX 	uintptr_t 	0 	UINTPTR_MAX

    return 0; 
} 
