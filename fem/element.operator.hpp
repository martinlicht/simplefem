#ifndef INCLUDEGUARD_FEM_ELEMENT_MATRIX
#define INCLUDEGUARD_FEM_ELEMENT_MATRIX


#include <vector>
#include <iostream>
// #include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/densematrix.hpp"
#include "../operators/linearoperator.hpp"







/* 
 * 
 * Transfer between matrix index and logical index *
 * 
 * /


/* calP_r\Lambda^k */

int logicalindex_into_matrixindex( int n, int r, int k, int& full_i, int  poly_i, int  sigma_i );

int logicalindex_from_matrixindex( int n, int r, int k, int  full_i, int& poly_i, int& sigma_i );

/* \mathring calP_r\Lambda^k */

int mathring_logicalindex_into_matrixindex( int n, int r, int k, int& full_i, int  poly_i, int  sigma_i );

int mathring_logicalindex_from_matrixindex( int n, int r, int k, int  full_i, int& poly_i, int& sigma_i );

/* calP_r^-\Lambda^k */

int trimmed_logicalindex_into_matrixindex( int n, int r, int k, int& full_i, int  poly_i, int  sigma_i );

int trimmed_logicalindex_from_matrixindex( int n, int r, int k, int  full_i, int& poly_i, int& sigma_i );

/* \mathring calP_r^-\Lambda^k */

int mathring_trimmed_logicalindex_into_matrixindex( int n, int r, int k, int& full_i, int  poly_i, int  sigma_i );

int mathring_trimmed_logicalindex_from_matrixindex( int n, int r, int k, int  full_i, int& poly_i, int& sigma_i );



/* Basis calP_r\Lambda^k */

int basis_logicalindex_into_matrixindex( int n, int r, int k, int& full_i, int  poly_i, int  sigma_i );

int basis_logicalindex_from_matrixindex( int n, int r, int k, int  full_i, int& poly_i, int& sigma_i );

/* Basis \mathring calP_r\Lambda^k */

int basis_mathring_logicalindex_into_matrixindex( int n, int r, int k, int& full_i, int  poly_i, int  sigma_i );

int basis_mathring_logicalindex_from_matrixindex( int n, int r, int k, int  full_i, int& poly_i, int& sigma_i );

/* Basis calP_r^-\Lambda^k */

int basis_trimmed_logicalindex_into_matrixindex( int n, int r, int k, int& full_i, int  poly_i, int  sigma_i );

int basis_trimmed_logicalindex_from_matrixindex( int n, int r, int k, int  full_i, int& poly_i, int& sigma_i );

/* Basis \mathring calP_r^-\Lambda^k */

int basis_mathring_trimmed_logicalindex_into_matrixindex( int n, int r, int k, int& full_i, int  poly_i, int  sigma_i );

int basis_mathring_trimmed_logicalindex_from_matrixindex( int n, int r, int k, int  full_i, int& poly_i, int& sigma_i );







/* 
 * 
 * Inclusions into larger spaces:
 * 
 * full poly -> next Whitney
 * Whitney -> next full_poly
 * mathring -> no-mathring
 * 
 * basis -> non-basis
 * 
 */

/* Degree elevation : steps as small as possible */

DenseMatrix inclusionmatrix_elevation_Polydiff_to_Whitney(
                int innerdimension, 
                int polynomialdegree, 
                int form degree 
                );

DenseMatrix inclusionmatrix_elevation_Whitney_to_Polydiff(
                int innerdimension, 
                int polynomialdegree, 
                int form degree 
                );

DenseMatrix inclusionmatrix_elevation_mathring_Polydiff_to_Whitney(
                int innerdimension, 
                int polynomialdegree, 
                int form degree 
                );

DenseMatrix inclusionmatrix_elevation_mathring_Whitney_to_Polydiff(
                int innerdimension, 
                int polynomialdegree, 
                int form degree 
                );


/* Skip the boundary conditions */

DenseMatrix inclusionmatrix_mathring_Polydiff(
                int innerdimension, 
                int polynomialdegree, 
                int form degree 
                );

DenseMatrix inclusionmatrix_mathring_Whitney(
                int innerdimension, 
                int polynomialdegree, 
                int form degree 
                );


/* From Basis to no Basis */

DenseMatrix inclusionmatrix_basisnobasis_Polydiff(
                int innerdimension, 
                int polynomialdegree, 
                int form degree 
                );

DenseMatrix inclusionmatrix_basisnobasis_Whitney(
                int innerdimension, 
                int polynomialdegree, 
                int form degree 
                );

DenseMatrix inclusionmatrix_basisnobasis_mathring_Polydiff(
                int innerdimension, 
                int polynomialdegree, 
                int form degree 
                );

DenseMatrix inclusionmatrix_basisnobasis_mathring_Whitney(
                int innerdimension, 
                int polynomialdegree, 
                int form degree 
                );





/* 
 * 
 * Trace Matrix 
 * 
 */

DenseMatrix element_trace(
                int innerdimension, 
                int polynomialdegree, 
                int form degree,
                int subsimplex_dim,
                int subsim_index
                );


/* 
 * 
 * Extension Matrix 
 * - uses the special extension operators
 *   discussed by Arnold, Falk, and Winther
 * 
 */

DenseMatrix element_extension(
                int innerdimension, 
                int polynomialdegree, 
                int form degree,
                int subsimplex_dim,
                int subsim_index
                );


/* 
 * 
 * Differential Matrix 
 * 
 */

DenseMatrix element_exteriorderivative(
                int innerdimension, 
                int polynomialdegree, 
                int form degree 
                );


/* 
 * 
 * Mass Matrix 
 * 
 */

DenseMatrix element_massproduct(
                int innerdimension, 
                int polynomialdegree, 
                int form degree,
                DenseMatrix distances
                );









#endif