
#include <ostream>
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

#include "../fem/global.unphysical.hpp"


SparseMatrix FEECCanonicalizeBroken( const Mesh& mesh, int n, int k, int r )
{
    Assert( 0 <= n );
    Assert( n <= mesh.getinnerdimension() );
    Assert( 0 <= r );
    Assert( 0 <= k and k <= n );


    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int poly_size = binomial_integer( n+r, n );
    const int form_size = binomial_integer( n+1, k );

    const int localdim = poly_size * form_size;
    
    const int dim_in      = num_simplices * localdim;
    const int dim_out     = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    

    // Calculate local matrix 

    DenseMatrix Aux1( n+1, n+1, 0. );
    for( int i = 1; i <= n; i++ ) {
        Aux1(i,i) = 1.;
        Aux1(i,0) = -1.;
    }
    
    const DenseMatrix Aux2 = SubdeterminantMatrix( Aux1, k );

    assert( Aux2.issquare() and Aux2.getdimout() == form_size );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    for( int i = 0; i < localdim; i++ )
    for( int j = 0; j < localdim; j++ )
    {
        int index_of_entry = s * localdim * localdim + i * localdim + j;

        int p_row = i / form_size;
        int p_col = j / form_size;

        int f_row = i % form_size;
        int f_col = j % form_size;
        
        SparseMatrix::MatrixEntry entry;
        entry.row    = s * localdim + i;
        entry.column = s * localdim + j;
        entry.value  = ( p_row == p_col ) ? Aux2( f_row, f_col ) : 0.;
        
        ret.setentry( index_of_entry, entry );
    }

    
    return ret;

}


SparseMatrix FEECRandomizeBroken( const Mesh& mesh, int n, int k, int r )
{
    Assert( 0 <= n );
    Assert( n <= mesh.getinnerdimension() );
    Assert( 0 <= r );
    Assert( 0 <= k and k <= n );


    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int poly_size = binomial_integer( n+r, n );
    const int form_size = binomial_integer( n+1, k );

    const int localdim = poly_size * form_size;
    
    const int dim_in      = num_simplices * localdim;
    const int dim_out     = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    

    // Calculate local matrix 

    int t = random_integer() % (n+1);
    
    // int alpha = random_uniform();
    Float alpha = random_uniform(); //static_cast<Float>( rand() ) / static_cast<Float>( RAND_MAX );
    Assert( 0. <= alpha and alpha <= 1. );

    DenseMatrix Aux1( n+1, n+1, 0. );
    for( int i = 0; i < t; i++ ) {
        Aux1(i,i) =  1.;
        Aux1(i,t) = -alpha;
    }
    Aux1(t,t) = 1. - alpha;    
    for( int i = t+1; i <= n; i++ ) {
        Aux1(i,i) = 1.;
        Aux1(i,t) = -alpha;
    }
    
    const DenseMatrix Aux2 = SubdeterminantMatrix( Aux1, k );

    assert( Aux2.issquare() and Aux2.getdimout() == form_size );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    for( int i = 0; i < localdim; i++ )
    for( int j = 0; j < localdim; j++ )
    {
        int index_of_entry = s * localdim * localdim + i * localdim + j;

        int p_row = i / form_size;
        int p_col = j / form_size;

        int f_row = i % form_size;
        int f_col = j % form_size;
        
        SparseMatrix::MatrixEntry entry;
        entry.row    = s * localdim + i;
        entry.column = s * localdim + j;
        entry.value  = ( p_row == p_col ) ? Aux2( f_row, f_col ) : 0.;
        
        ret.setentry( index_of_entry, entry );
    }

    
    return ret;

}
