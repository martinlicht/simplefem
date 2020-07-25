#ifndef INCLUDEGUARD_FEM_FEECBROKENELEVATIONMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENELEVATIONMATRIX


#include <vector>
#include <iostream>
// #include <cassert>

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
//  Matrix for degree elevantion                    //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//  and lifts up rplus degrees                      //
//                                                  //
//////////////////////////////////////////////////////





inline SparseMatrix FEECBrokenElevationMatrix( Mesh& mesh, int n, int k, int r, int rplus = 1 )
{
    
    // check whether the parameters are right 
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( rplus >= 1 );
    
    
    
    const std::vector<MultiIndex> multis_adds = generateMultiIndices( IndexRange( 0, n ), rplus     );
    const std::vector<MultiIndex> multis_low  = generateMultiIndices( IndexRange( 0, n ), r         );
    const std::vector<MultiIndex> multis_high = generateMultiIndices( IndexRange( 0, n ), r + rplus );
    
    assert( multis_adds.size() == binomial_integer( n + rplus    , n ) );
    assert( multis_low.size()  == binomial_integer( n + r,         n ) );
    assert( multis_high.size() == binomial_integer( n + r + rplus, n ) );
    
    const std::vector<IndexMap> sigmas = generateSigmas( IndexRange( 1, k ), IndexRange( 0, n ) );
    
    assert( sigmas.size() == binomial_integer( n+1, k ) );
    
    const int localdim_in  = multis_low.size()  * sigmas.size();
    const int localdim_out = multis_high.size() * sigmas.size();
    
    std::vector<SparseMatrix::MatrixEntry> localmatrixentries;
    
    for( int low_poly_index = 0; low_poly_index < multis_low.size(); low_poly_index++ )
    for( int     form_index = 0;     form_index <     sigmas.size();     form_index++ )
    for( const MultiIndex& addendum : multis_adds )
    {
        
        const MultiIndex& low_poly = multis_low[low_poly_index];
        
        MultiIndex high_poly = low_poly + addendum;
        
        int high_poly_index = find_index( multis_high, high_poly );
        
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = high_poly_index * sigmas.size() + form_index;
        entry.column = low_poly_index * sigmas.size() + form_index;
        entry.value  = 1.0;
        
        assert( entry.row >= 0 && entry.row < localdim_out );
        assert( entry.column >= 0 && entry.column < localdim_in );
        
        localmatrixentries.push_back( entry );
        
    }
    
    // Auxiliary calculations and preparations
    
    // Finished generating local matrix
    
    const int num_simplices = mesh.count_simplices( n );
        
    int noe = localmatrixentries.size();
    
    const int dim_out = num_simplices * localdim_out;
    const int dim_in  = num_simplices * localdim_in;
    
    const int num_entries = num_simplices * noe;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices;   s++ )
    for( int i = 0; i < noe; i++ )
    {
        
        int index_of_entry = s * noe + i; 
            
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = s * localdim_out + localmatrixentries[i].row;
        entry.column = s * localdim_in  + localmatrixentries[i].column;
        entry.value  = localmatrixentries[i].value;
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}




#endif
