
#include <vector>
#include <iostream>
// #include <cassert>

#include "../basic.hpp"
#include "../utility/stl.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/global.trace.hpp"




SparseMatrix FEECBrokenTraceMatrix( const Mesh& mesh, int n, int k, int r, bool is_signed )
{
    
    Assert( 0 <= n );
    Assert( n <= mesh.getinnerdimension() );
    Assert( 0 <= r );
    Assert( 0 <= k and k <= n );

    // TODO: catch the border cases r <= 0 or k >= n
    
    const std::vector<MultiIndex>& cell_multis = generateMultiIndices( IndexRange( 0, n   ), r );
    const std::vector<MultiIndex>& face_multis = generateMultiIndices( IndexRange( 0, n-1 ), r );
    
    const int dim_cell_polynomials = cell_multis.size();
    const int dim_face_polynomials = face_multis.size();

    const std::vector<IndexMap>& cell_sigmas = generateSigmas( IndexRange( 1, k ), IndexRange( 0, n   ) );
    const std::vector<IndexMap>& face_sigmas = generateSigmas( IndexRange( 1, k ), IndexRange( 0, n-1 ) );
    
    const int dim_cell_sigmas = cell_sigmas.size();
    const int dim_face_sigmas = face_sigmas.size();
    
    const int num_cells = mesh.count_simplices(n  );
    const int num_faces = mesh.count_simplices(n-1);

    const int dim_in      = num_cells * dim_cell_polynomials * dim_cell_sigmas;
    const int dim_out     = num_faces * dim_face_polynomials * dim_face_sigmas;
    
    const int num_entries = num_cells * (n+1) * dim_face_polynomials * dim_face_sigmas;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );

    auto face_inclusions = generateSigmas( IndexRange(0,n-1), IndexRange(0,n) );
    

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s  = 0; s  <  num_cells;  s++ )
    for( int fi = 0; fi <=         n; fi++ )
    {
        
        // Find the inclusion of the face f in the cell s
        
        const auto vertices_of_cell = mesh.getsubsimplices( n, 0, s );

        const int face = mesh.get_subsimplex( n, n-1, s, fi );
        
        const auto vertices_of_face = mesh.getsubsimplices( n-1, 0, face );

        // assert( vertices_of_cell.isstrictlyascending() );
        // assert( vertices_of_face.isstrictlyascending() );

        int inclusion_index = 0;
        for( ; inclusion_index < face_inclusions.size(); inclusion_index++ )
            if( vertices_of_face == vertices_of_cell * face_inclusions[inclusion_index] )
                break;
        assert( inclusion_index < face_inclusions.size() );
        assert( vertices_of_face == vertices_of_cell * face_inclusions[inclusion_index] );
        const IndexMap inclusion_face_to_cell = face_inclusions[inclusion_index];
        

        int index_of_gap = 0;
        for( int g = 0; g <= n-1; g++ )
            if( inclusion_face_to_cell[g] != g ) {
                index_of_gap = g;
                break;
            }
        const int signum = signpower( index_of_gap ); // TODO: this is incorrect


        for( int mi = 0; mi < face_multis.size(); mi++ )
        for( int fs = 0; fs < face_sigmas.size(); fs++ )
        {

            const MultiIndex& face_multi = face_multis[mi];

            const auto face_multi_index = mi;

            const IndexMap& face_sigma = face_sigmas[fs];

            const auto face_sigma_index = fs;

            const MultiIndex cell_multi = MultiIndex( IndexRange(0,n), [&face_multi,&inclusion_face_to_cell]( int p ) -> int {
                                            assert( inclusion_face_to_cell.getTargetRange().contains(p) ); 
                                            if( inclusion_face_to_cell.rangecontains(p) )
                                                return face_multi.at( inclusion_face_to_cell.preimageof(p) );
                                            else
                                                return 0;
                                        } );
            
            int cell_multi_index = find_index( cell_multis, cell_multi );
                
            const IndexMap cell_sigma = inclusion_face_to_cell * face_sigma;

            int cell_sigma_index = find_index( cell_sigmas, cell_sigma );


            // Having found the relevant local indices, we can now construct the global indices

            SparseMatrix::MatrixEntry entry;
        
            entry.row    = face * dim_face_polynomials * dim_face_sigmas + face_multi_index * dim_face_sigmas + face_sigma_index;
            entry.column = s    * dim_cell_polynomials * dim_cell_sigmas + cell_multi_index * dim_cell_sigmas + cell_sigma_index;
            entry.value  = is_signed ? signum : 1.0;
            

            int index_of_entry = s * ( (n+1) * dim_face_polynomials * dim_face_sigmas ) + fi * ( dim_face_polynomials * dim_face_sigmas ) + face_multi_index * dim_face_sigmas + face_sigma_index;
            // num_cells * (n+1) * dim_face_polynomials * dim_face_sigmas;

            ret.setentry( index_of_entry, entry );
            
        }
    
    }

    assert( ret.isfinite() );
    return ret;
    
}



