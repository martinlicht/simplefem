#ifndef INCLUDEGUARD_SOLVER_SYSTEMSOLVER
#define INCLUDEGUARD_SOLVER_SYSTEMSOLVER

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <new>
#include <utility>

#include "../basic.hpp"
#include "sparsesolver.hpp"

const Float expected_sign_of_A =  1.;

void HodgeConjugateResidualSolverCSR( 
    const int N, 
    const int L, 
    Float* x, 
    const Float* b, 
    const int*  Arows, const int*  Acolumns, const Float*  Avalues, 
    const int*  Brows, const int*  Bcolumns, const Float*  Bvalues, 
    const int* Btrows, const int* Btcolumns, const Float* Btvalues, 
    const int*  Crows, const int*  Ccolumns, const Float*  Cvalues, 
    Float* residual,
    Float threshold,
    int print_modulo,
    Float inneriteration_threshold,
    int inneriteration_print_modulo
);

void HodgeConjugateResidualSolverCSR_diagonal( 
    const int N, 
    const int L, 
    Float* x, 
    const Float* b, 
    const int*  Arows, const int*  Acolumns, const Float*  Avalues, 
    const int*  Brows, const int*  Bcolumns, const Float*  Bvalues, 
    const int* Btrows, const int* Btcolumns, const Float* Btvalues, 
    const int*  Crows, const int*  Ccolumns, const Float*  Cvalues, 
    Float* residual,
    Float threshold,
    int print_modulo,
    Float inneriteration_threshold,
    int inneriteration_print_modulo
);

void HodgeConjugateResidualSolverCSR_textbook( 
    const int N, 
    const int L, 
    Float* x, 
    const Float* b, 
    const int*  Arows, const int*  Acolumns, const Float*  Avalues, 
    const int*  Brows, const int*  Bcolumns, const Float*  Bvalues, 
    const int* Btrows, const int* Btcolumns, const Float* Btvalues, 
    const int*  Crows, const int*  Ccolumns, const Float*  Cvalues, 
    Float* residual,
    Float threshold,
    int print_modulo,
    Float inneriteration_threshold,
    int inneriteration_print_modulo
);

void HodgeConjugateResidualSolverCSR_SSOR( 
    const int N, 
    const int L, 
    Float* x, 
    const Float* b, 
    const int*  Arows, const int*  Acolumns, const Float*  Avalues, 
    const int*  Brows, const int*  Bcolumns, const Float*  Bvalues, 
    const int* Btrows, const int* Btcolumns, const Float* Btvalues, 
    const int*  Crows, const int*  Ccolumns, const Float*  Cvalues, 
    Float* residual,
    Float threshold,
    int print_modulo,
    Float inneriteration_threshold,
    int inneriteration_print_modulo
);


 
#endif
