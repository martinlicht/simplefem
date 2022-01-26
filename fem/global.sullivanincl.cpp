
#include <iostream>
#include <utility>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/linearoperator.hpp"
#include "../mesh/mesh.hpp"

#include "indexfunctions.hpp"

#include "../fem/global.sullivanincl.hpp"


/* 
TODO:
- The major difficulty here is whether the subsimplex indices of the vertices are ordered or not.
- The lack of a general ordering principle forecloses an efficient implmentation
where the loop is completely streamed.
- 
*/

SparseMatrix FEECSullivanInclusionMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( r >= 1 or k == n ); // TODO: wir erlauben r==0 im falle von volumenformen 

    // generate the list of sigmas and multiindices 
    
    std::vector< std::vector< std::pair<MultiIndex,IndexMap> > > lists_of_sullivan_indices( n+1 );
    for( auto dimension : IndexRange(0,n) )
        lists_of_sullivan_indices[dimension] = ListOfSullivanIndices( dimension, k, r );

    const auto multis_dst = generateMultiIndices( IndexRange(0,n), r );
    const auto sigmas_dst = generateSigmas( IndexRange(1,k), IndexRange(0,n) );
    
    
    // count mesh properties 
    
    const int num_volumes = mesh.count_simplices( n );
    
    std::vector< int > num_faces( n+1 );
    for( auto dimension : IndexRange(0,n) )
        num_faces[dimension] = mesh.count_simplices( dimension );

    
    // dimensions of each global matrix 
    
    const int dim_out = num_volumes * binomial_integer( n+1, k ) * binomial_integer( n+r, r );
    
    // 2+1+1 choose 1
    // const int dim_in = sum_int( n, [ num_faces&, lists_of_sullivan_indices& ](int d) 
    //                      -> int{ return num_faces[d] * lists_of_sullivan_indices[d].size(); } )
    
    int dim_in  = 0;
    for( auto d : IndexRange(0,n) )
        dim_in += num_faces[d] * lists_of_sullivan_indices[d].size();
    
    // const int num_entries = num_volumes * sum_int( n, [ num_faces&, lists_of_sullivan_indices&, n ](int d) -> int{ return binomial_integer(n+1,d+1) * lists_of_sullivan_indices[d].size(); } )
    int num_entries = 0;
    for( auto d : IndexRange(0,n) )
        num_entries += num_volumes * binomial_integer(n+1,d+1) * lists_of_sullivan_indices[d].size();
    
    // create that sparse matrix 
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    
    // for auxiliary purposes, create the index inclusion maps 
    // from subsimplices into the supersimplices
    
    std::vector< std::vector<IndexMap> > subsimplex_inclusions( n+1 );
    for( auto d : IndexRange(0,n) )
        subsimplex_inclusions[d] = generateSigmas( IndexRange(0,d), IndexRange(0,n) );
    
    // go over all the subsimplex dimensions
    // go over all the d dimensional reference subsimplices 
    // go over the corresponding alpha/sigma pairs
    // go over all the volumes 

    for( int d  = 0; d  <= n;                          d++ )       
    for( int fi = 0; fi  < binomial_integer(n+1,d+1); fi++ )       
    // for( int index_alphasigma = 0; index_alphasigma < lists_of_sullivan_indices[d]; index_alphasigma++ )
    for( const auto& alphasigma : lists_of_sullivan_indices[d] )   
    for( int s  = 0; s   < num_volumes;                s++ )       
    {
        
        // Find indices of things and prepare auxiliary variables 
        
        const MultiIndex& alpha = alphasigma.first;
        
        const IndexMap&   sigma = alphasigma.second;

        // TODO: This search can be replaced by just changing the loop above for f's sake 
        // use the commented out for loop above 

        const int index_alphasigma = 
            std::find( lists_of_sullivan_indices[d].begin(), lists_of_sullivan_indices[d].end(), alphasigma )
            - 
            lists_of_sullivan_indices[d].begin();
        assert( 0 <= index_alphasigma && index_alphasigma < lists_of_sullivan_indices[d].size() );
        assert( lists_of_sullivan_indices[d][index_alphasigma] == alphasigma );
        
        const int index_fi = mesh.get_subsimplex( n, d, s, fi );
        assert( 0 <= index_fi && index_fi < mesh.count_simplices(d) );
        
        
        // get inclusion index map
                
        IndexMap volume_vertices = mesh.getsubsimplices( n, 0, s        );
        IndexMap face_vertices   = mesh.getsubsimplices( d, 0, index_fi );
        int inclusion_index = 0;
        for( ; inclusion_index < subsimplex_inclusions[d].size(); inclusion_index++ )
            if( face_vertices == volume_vertices * subsimplex_inclusions[d][inclusion_index] )
                break;
        assert( inclusion_index < subsimplex_inclusions[d].size() );
        assert( face_vertices == volume_vertices * subsimplex_inclusions[d][inclusion_index] );
        const IndexMap inclusion = subsimplex_inclusions[d][inclusion_index];
        
        
        // create actual multiindices 
        
        const MultiIndex alpha_vol = MultiIndex( IndexRange(0,n), [&alpha,&inclusion]( int p ) -> int {
                                            assert( inclusion.getTargetRange().contains(p) ); 
                                            if( inclusion.rangecontains(p) )
                                                return alpha.at( inclusion.preimageof(p) );
                                            else
                                                return 0;
                                        } );
        
        const IndexMap   sigma_vol = inclusion * sigma;
        
        
        // find the indices of those extended functions 
        
        const auto index_alpha_vol = std::find( multis_dst.begin(), multis_dst.end(), alpha_vol ) - multis_dst.begin();
        
        const auto index_sigma_vol = std::find( sigmas_dst.begin(), sigmas_dst.end(), sigma_vol ) - sigmas_dst.begin();

        assert( 0 <= index_alpha_vol and index_alpha_vol < multis_dst.size() );
        assert( 0 <= index_sigma_vol and index_sigma_vol < sigmas_dst.size() );
        

        // enter the values of the data structure 

        int rowindex = s * binomial_integer(n+r,r) * binomial_integer(n+1,k) 
                       +
                       index_alpha_vol * binomial_integer(n+1,k)
                       + 
                       index_sigma_vol;

        int colindex = sum_int( d-1, 
                        [&mesh,&lists_of_sullivan_indices](int i) -> int { return mesh.count_simplices(i) * lists_of_sullivan_indices[i].size(); } 
                       )
                       + 
                       index_fi * lists_of_sullivan_indices[d].size()
                       +
                       index_alphasigma
                       ;

        Float value  = 1.;

        int index_of_entry = sum_int( d-1, 
                                [ &lists_of_sullivan_indices, n ](int c) -> int { return binomial_integer(n+1,c+1) * lists_of_sullivan_indices[c].size(); }
                             ) * num_volumes
                             +
                             fi * lists_of_sullivan_indices[d].size() * num_volumes
                             +
                             index_alphasigma * num_volumes
                             +
                             s; 
                             
        assert( rowindex       < dim_out );
        assert( colindex       < dim_in  );
        assert( index_of_entry < num_entries );
        assert( rowindex >= 0 && colindex >= 0 && index_of_entry >= 0 );
            
        // set up the actual entry
        
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = rowindex;
        entry.column = colindex;
        entry.value  = value;
        
        if( mesh.get_flag( d, index_fi ) == SimplexFlagDirichlet )
            entry.value = 0.;
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}

