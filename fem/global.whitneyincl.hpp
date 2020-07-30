#ifndef INCLUDEGUARD_FEM_FEECWHITNEYINCLUSIONMATRIX
#define INCLUDEGUARD_FEM_FEECWHITNEYINCLUSIONMATRIX


// #include <cassert>
#include <iostream>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/matrixtensorproduct.hpp"
#include "../operators/linearoperator.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/local.polynomialmassmatrix.hpp"



//////////////////////////////////////////////////////
//                                                  //
//  Matrix for inclusion of                         //
//  continuous whitney forms                        //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//  and lifts up rplus degrees                      //
//                                                  //
//////////////////////////////////////////////////////




inline SparseMatrix FEECWhitneyInclusionMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( r >= 1 );
    assert( r == 1 ); // only lowest order supported
    
    const auto sigmas_src = generateSigmas( IndexRange(0,k  ), IndexRange(0,n) );
    const auto sigmas_dst = generateSigmas( IndexRange(0,k-1), IndexRange(0,n) );
    
    const int num_volumes = mesh.count_simplices( n );
    const int num_faces   = mesh.count_simplices( k );
    
    const int dim_out = num_volumes * binomial_integer( n+1, k ) * (n+1);
    const int dim_in  = num_faces;
    
    const int num_entries = num_volumes * binomial_integer( n+1, k+1 ) * (k+1);
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s  = 0; s  < num_volumes;                  s++  )
    for( int fi = 0; fi < binomial_integer( n+1, k+1 ); fi++ )
    for( int p  = 0; p <= k                           ; p++  )
    {
        
        int sigma_src_index = -1; // TODO find the sigma corresponding to fi
        IndexMap sigma_src = sigmas_src[ sigma_src_index ];

        IndexMap sigma_dst = sigma_src - sigma_src[p]; // TODO: pop the p-th element 
        int sigma_dst_index = find_index( sigmas_dst, sigma_dst ); // TODO the new sigma in the 
        
        int rowindex = p * binomial_integer( n+1, k ) + sigma_dst_index;

        int colindex = mesh.get_subsimplex( n, k, s, fi );

        Float value  = signpower(p);

        int index_of_entry = s * binomial_integer( n+1, k+1 ) * (k+1) + fi * binomial_integer( n+1, k+1 ) + p; 
            
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = rowindex;
        entry.column = colindex;
        entry.value  = signpower(p);
        
        if( mesh.get_flag( k, mesh.get_subsimplex( n, k, s, fi ) ) == SimplexFlagDirichlet )
            entry.value = 0.;
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}




#endif
