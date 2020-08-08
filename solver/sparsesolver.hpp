#ifndef INCLUDEGUARD_SOLVER_CLASSICALSOLVERS
#define INCLUDEGUARD_SOLVER_CLASSICALSOLVERS


#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <memory.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include "../basic.hpp"



// Solves the sparse matrix system with conjugate gradients

void ConjugateGradientSolverCSR( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float allowed_error,
    unsigned int restart_modulo
);


// Solves the sparse matrix system with conjugate residuals

void ConjugateResidualSolverCSR( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float allowed_error,
    unsigned int restart_modulo
);



// // Solves the Uzawa-system
// 
// void UzawaConjugateResidualCSR(
//     const int M, const int N,
//     Float* const sigma, Float* const u,
//     Float* const ressigma, Float *const resu,
//     const Float* const e, const Float* const f,
//     const Float* const entriesA, const int* const csrrowsA, const int* const csrcolumnsA,
//     const Float* const entriesB, const int* const csrrowsB, const int* const csrcolumnsB,
//     const Float* const entriesC, const int* const csrrowsC, const int* const csrcolumnsC
// );
                            

#endif
