
#include <utility>
#include <vector>

#include "../basic.hpp"
#include "../utility/stl.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/linearoperator.hpp"
#include "../mesh/mesh.hpp"

#include "indexfunctions.hpp"

#include "../fem/global.avgsullivan.hpp"


/* 
TODO:
- The major difficulty here is whether the subsimplex indices of the vertices are ordered or not.
- The lack of a general ordering principle forecloses an efficient implmentation where the loop is completely uniform.
*/

SparseMatrix FEECSullivanAveragingMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( r >= 1 or k == n ); // TODO: we allow r == 0 in case of volume forms 

    // generate the list of sigmas and multiindices for each dimension 
    
    std::vector< std::vector< std::pair<MultiIndex,IndexMap> > > lists_of_sullivan_indices( n+1 );
    for( auto d : IndexRange(0,n) )
        lists_of_sullivan_indices[d] = ListOfSullivanIndices( d, k, r );

    const auto multis_dst = generateMultiIndices( IndexRange(0,n), r );
    const auto sigmas_dst = generateSigmas( IndexRange(1,k), IndexRange(0,n) );
    
    
    // count mesh properties 
    
    const int num_volumes = mesh.count_simplices( n );
    
    std::vector< int > num_faces( n+1 );
    for( auto d : IndexRange(0,n) )
        num_faces[d] = mesh.count_simplices( d );

    
    // dimensions of the matrix 
    
    const int dim_in = num_volumes * binomial_integer( n+1, k ) * binomial_integer( n+r, r );
    
    int dim_out  = 0;
    for( auto d : IndexRange(0,n) )
        dim_out += num_faces[d] * lists_of_sullivan_indices[d].size();
    
    int num_entries = 0;
    for( auto d : IndexRange(0,n) )
        num_entries += num_volumes * binomial_integer(n+1,d+1) * lists_of_sullivan_indices[d].size();
    
    // create that matrix 
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    
    // for auxiliary purposes, create the index inclusion maps 
    // from subsimplices into the supersimplices
    
    std::vector< std::vector<IndexMap> > subsimplex_inclusions( n+1 );
    for( auto d : IndexRange(0,n) )
        subsimplex_inclusions[d] = generateSigmas( IndexRange(0,d), IndexRange(0,n) );
    
    for( int d  = 0; d  <= n;                          d++ )                                                      // go over all the subsimplex dimensions
    for( int fi = 0; fi  < binomial_integer(n+1,d+1); fi++ )                                                      // go over all the d dimensional subsimplices 
    for( int index_alphasigma = 0; index_alphasigma < lists_of_sullivan_indices[d].size(); index_alphasigma++ )   // go over the corresponding alpha/sigma pairs
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif 
    for( int s  = 0; s   < num_volumes;                s++ )                                                      // go over all the volumes 
    {
        
        // Find indices of things and prepare auxiliary variables 
        
        const auto& alphasigma = lists_of_sullivan_indices[d][index_alphasigma];

        const MultiIndex& alpha = alphasigma.first;
        
        const IndexMap&   sigma = alphasigma.second;
        
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

        int colindex = s * binomial_integer(n+r,r) * binomial_integer(n+1,k) 
                       +
                       index_alpha_vol * binomial_integer(n+1,k)
                       + 
                       index_sigma_vol;

        int rowindex = sum_int( d-1, 
                        [&mesh,&lists_of_sullivan_indices](int i) -> int { return mesh.count_simplices(i) * lists_of_sullivan_indices[i].size(); } 
                       )
                       + 
                       index_fi * lists_of_sullivan_indices[d].size()
                       +
                       index_alphasigma
                       ;

        Float value  = 1. / mesh.getsupersimplices( n, d, index_fi ).size();

        int index_of_entry = sum_int( d-1, 
                                [ &lists_of_sullivan_indices, n ](int c) -> int { return binomial_integer(n+1,c+1) * lists_of_sullivan_indices[c].size(); }
                             ) * num_volumes
                             +
                             fi * lists_of_sullivan_indices[d].size() * num_volumes
                             +
                             index_alphasigma * num_volumes
                             +
                             s; 
                             
        assert( rowindex >= 0 && colindex >= 0 && index_of_entry >= 0 );
        assert( rowindex       < dim_out );
        assert( colindex       < dim_in  );
        assert( index_of_entry < num_entries );
            
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

