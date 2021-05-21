#ifndef INCLUDEGUARD_SOLVER_CLASSICALSOLVERS
#define INCLUDEGUARD_SOLVER_CLASSICALSOLVERS


#include <new>

#include "../basic.hpp"



// Solves the sparse matrix system with conjugate gradients

void ConjugateGradientSolverCSR( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float threshold,
    int print_modulo
);

void ConjugateGradientSolverCSR_DiagonalPreconditioner( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float threshold,
    int print_modulo,
    const Float* precon 
);


void ConjugateGradientSolverCSR_SSOR( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float threshold,
    int print_modulo,
    const Float* precon,
    Float omega
);


// Solves the sparse matrix system with conjugate residuals

void ConjugateResidualSolverCSR( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float threshold,
    int print_modulo
);

void ConjugateResidualSolverCSR_textbook( 
    const int N, 
    Float* x, 
    const Float* b, 
    const int* csrrows, const int* csrcolumns, const Float* csrvalues, 
    Float* residual,
    Float threshold,
    int print_modulo
);




void MINRESCSR( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ res,
    const Float threshold,
    int print_modulo
);



void WHATEVER( 
    const int N, 
    Float* __restrict__ x, 
    const Float* __restrict__ b, 
    const int* __restrict__ csrrows, const int* __restrict__ csrcolumns, const Float* __restrict__ csrvalues, 
    Float* __restrict__ res,
    const Float threshold,
    int print_modulo
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
