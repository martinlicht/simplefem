#include <vector>

#include "../base/include.hpp"
#include "../utility/stl.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"
#include "../utility/random.hpp"

#include "global.diffmatrix.hpp"
#include "global.elevation.hpp"
#include "global.interpol.hpp"
#include "global.unphysical.hpp"
#include "utilities.hpp"






SparseMatrix FEECBrokenDiffMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 1 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <  n );
    
    
    // generate an empty local sample matrix
    
    const std::vector<MultiIndex> multis_dest = generateMultiIndices( IndexRange( 0, n ), r-1 );
    const std::vector<MultiIndex> multis_src  = generateMultiIndices( IndexRange( 0, n ), r );
    
    const std::vector<IndexMap>   sigmas_dest = generateSigmas( IndexRange( 1, k+1 ), IndexRange( 0, n ) );
    const std::vector<IndexMap>   sigmas_src  = generateSigmas( IndexRange( 1, k   ), IndexRange( 0, n ) );
    
    const int localdim_out = multis_dest.size() * sigmas_dest.size();
    const int localdim_in  = multis_src.size()  * sigmas_src.size();
    
    std::vector<SparseMatrix::MatrixEntry> localmatrixentries;
    
    assert( multis_dest.size() == binomial_integer( n+r-1, n   ) );
    assert( multis_src.size()  == binomial_integer( n+r,   n   ) );
    assert( sigmas_dest.size() == binomial_integer( n+1,   k+1 ) );
    assert( sigmas_src.size()  == binomial_integer( n+1,   k   ) );
    
    for( int src_poly_index = 0; src_poly_index < multis_src.size(); src_poly_index++ )
    for( int src_form_index = 0; src_form_index < sigmas_src.size(); src_form_index++ )
    for( int p = 0; p <= n; p++ )
    {
        
        const MultiIndex& src_poly = multis_src[src_poly_index];
        const IndexMap&   src_form = sigmas_src[src_form_index];
        
        if( src_poly[p] == 0 or src_form.has_value_in_range(p) )
            continue;
        
        MultiIndex new_poly = src_poly - p;
        IndexMap   new_form = expand_one( src_form, p );
        
        int new_poly_index = find_index( multis_dest, new_poly );
        int new_form_index = find_index( sigmas_dest, new_form );
        
        assert( new_form.getSourceRange().min() == 1 );
        int signum = sign_power( new_form.get_preimage_of(p) - 1 );
        
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = new_poly_index * sigmas_dest.size() + new_form_index;
        entry.column = src_poly_index * sigmas_src.size()  + src_form_index;
        entry.value  = src_poly[p] * signum;
        
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

    const int dim_in      = num_simplices * poly_size * form_size;
    const int dim_out     = num_simplices * poly_size * form_size;
    const int num_entries = num_simplices * poly_size * form_size * form_size;

    SparseMatrix ret( dim_out, dim_in, num_entries );


    // Calculate local matrix

    DenseMatrix Aux1( n+1, n+1, 0. );
    for( int i = 1; i <= n; i++ ) {
        Aux1(i,i) = 1.;
        Aux1(i,0) = -1.;
    }
    // DenseMatrix Aux1 = IdentityMatrix(n+1) - DenseMatrix( n+1, n+1, 1./(n+1) );

    const DenseMatrix Aux = SubdeterminantMatrix( Aux1, k );

    assert( Aux.is_square() and Aux.getdimout() == form_size );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s  = 0;  s < num_simplices;  s++ )
    for( int p  = 0;  p <     poly_size;  p++ )
    for( int f1 = 0; f1 <     form_size; f1++ )
    for( int f2 = 0; f2 <     form_size; f2++ )
    {
        int index_of_entry = s * poly_size * form_size * form_size + p * form_size * form_size + f1 * form_size + f2;

        SparseMatrix::MatrixEntry entry;
        entry.row    = s * poly_size * form_size + p * form_size + f1;
        entry.column = s * poly_size * form_size + p * form_size + f2;
        entry.value  = Aux( f1, f2 );

        ret.setentry( index_of_entry, entry );
    }

    return ret;

}





SparseMatrix FEECRandomizeBroken( const Mesh& mesh, int n, int k, int r, Float base_alpha )
{
    Assert( 0 <= n );
    Assert( n <= mesh.getinnerdimension() );
    Assert( 0 <= r );
    Assert( 0 <= k and k <= n );


    // Auxiliary calculations and preparations

    const int num_simplices = mesh.count_simplices( n );

    const int poly_size = binomial_integer( n+r, n );
    const int form_size = binomial_integer( n+1, k );

    const int dim_in      = num_simplices * poly_size * form_size;
    const int dim_out     = num_simplices * poly_size * form_size;
    const int num_entries = num_simplices * poly_size * form_size * form_size;

    SparseMatrix ret( dim_out, dim_in, num_entries );


    std::vector<DenseMatrix> auxiliaries( n+1, DenseMatrix(form_size,form_size,notanumber) );

    // Calculate local matrix
    for( int t = 0; t <= n; t++ )
    {
        Float alpha = 1.;
        alpha = ( 0. <= base_alpha and base_alpha <= 1.0 ) ? base_alpha : random_uniform();
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

        assert( Aux2.is_square() and Aux2.getdimout() == form_size );

        auxiliaries[t] = Aux2;
    }

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s  = 0;  s < num_simplices;  s++ )
    {
        int t = random_integer() % (n+1);
        assert( 0 <= t and t <= n );

        for( int p  = 0;  p <     poly_size;  p++ )
        for( int f1 = 0; f1 <     form_size; f1++ )
        for( int f2 = 0; f2 <     form_size; f2++ )
        {
            int index_of_entry = s * poly_size * form_size * form_size + p * form_size * form_size + f1 * form_size + f2;


            SparseMatrix::MatrixEntry entry;
            entry.row    = s * poly_size * form_size + p * form_size + f1;
            entry.column = s * poly_size * form_size + p * form_size + f2;
            entry.value  = auxiliaries[t]( f1, f2 );

            ret.setentry( index_of_entry, entry );
        }
    }

    return ret;

}













































SparseMatrix FEECBrokenElevationMatrix( const Mesh& mesh, int n, int k, int r, int r_plus )
{

    // check whether the parameters are right

    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( r_plus >= 0 );

    const std::vector<MultiIndex> multis_add  = generateMultiIndices( IndexRange( 0, n ), r_plus     );
    const std::vector<MultiIndex> multis_low  = generateMultiIndices( IndexRange( 0, n ), r          );
    const std::vector<MultiIndex> multis_high = generateMultiIndices( IndexRange( 0, n ), r + r_plus );

    assert( multis_add.size()  == binomial_integer( n + r_plus    , n ) );
    assert( multis_low.size()  == binomial_integer( n + r,          n ) );
    assert( multis_high.size() == binomial_integer( n + r + r_plus, n ) );

    const std::vector<IndexMap> sigmas = generateSigmas( IndexRange( 1, k ), IndexRange( 0, n ) );

    assert( sigmas.size() == binomial_integer( n+1, k ) );

    const int localdim_in  = multis_low.size()  * sigmas.size();
    const int localdim_out = multis_high.size() * sigmas.size();

    std::vector<SparseMatrix::MatrixEntry> localmatrixentries;

    for( int low_poly_index = 0; low_poly_index < multis_low.size(); low_poly_index++ )
    for( int     form_index = 0;     form_index <     sigmas.size();     form_index++ )
    for( int add_poly_index = 0; add_poly_index < multis_add.size(); add_poly_index++ )
//     for( const MultiIndex& addendum : multis_add )
    {

        const MultiIndex& addendum = multis_add[add_poly_index];

        const MultiIndex& low_poly = multis_low[low_poly_index];

        MultiIndex high_poly = low_poly + addendum;

        int high_poly_index = find_index( multis_high, high_poly );

        assert( 0 <= high_poly_index and high_poly_index < multis_high.size() );

        SparseMatrix::MatrixEntry entry;

        entry.row    = high_poly_index * sigmas.size() + form_index;
        entry.column = low_poly_index  * sigmas.size() + form_index;

        entry.value  = factorial_numerical( r_plus ) / addendum.factorial_numerical();

        assert( entry.row    >= 0 && entry.row    < localdim_out );
        assert( entry.column >= 0 && entry.column < localdim_in  );

        localmatrixentries.push_back( entry );

    }

// //     if( k==0 and r==0 and r_plus == 3 ){
// //
// //         for( const auto& lme : localmatrixentries )
// //             LOG << lme.row << space << lme.column << space << lme.value;// << nl;
// //
// //         exit(0);
// //     }

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

























SparseMatrix FEECBrokenInterpolationMatrix( const Mesh& mesh, int n, int k, int r, int r_plus )
{
    
    // check whether the parameters are right 
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( r_plus >= 0 );
    
    const std::vector<MultiIndex> multis_low  = generateMultiIndices( IndexRange( 0, n ), r          );
    const std::vector<MultiIndex> multis_high = generateMultiIndices( IndexRange( 0, n ), r + r_plus );
    
    assert( multis_low.size()  == binomial_integer( n + r         , n ) );
    assert( multis_high.size() == binomial_integer( n + r + r_plus, n ) );
    
    const std::vector<IndexMap> sigmas = generateSigmas( IndexRange( 1, k ), IndexRange( 0, n ) );
    
    assert( sigmas.size() == binomial_integer( n+1, k ) );


    // check whether the parameters are right 

    const DenseMatrix bcs = InterpolationPointsInBarycentricCoordinates( n, r );

    const DenseMatrix pvom_plus = PointValuesOfMonomials( r + r_plus, bcs );

    const DenseMatrix pvom      = PointValuesOfMonomials( r         , bcs );

    assert( pvom.is_square() );

    const DenseMatrix lpc = LagrangePolynomialCoefficients( n, r );
    const DenseMatrix lpc_inv = Inverse(lpc);

    const DenseMatrix localmatrix = MatrixTensorProduct( Inverse( pvom ) * pvom_plus, IdentityMatrix( sigmas.size() ) );


    
    

    
    const int localdim_in  = multis_high.size() * sigmas.size();
    const int localdim_out = multis_low.size()  * sigmas.size();
    
    std::vector<SparseMatrix::MatrixEntry> localmatrixentries;
    
    for( int  low_poly_index = 0; low_poly_index  <  multis_low.size();  low_poly_index++ )
    for( int high_poly_index = 0; high_poly_index < multis_high.size(); high_poly_index++ )
    for( int      form_index = 0;      form_index <      sigmas.size();      form_index++ )
    {
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = low_poly_index  * sigmas.size() + form_index;
        entry.column = high_poly_index * sigmas.size() + form_index;
        entry.value  = localmatrix( entry.row, entry.column );
        
        assert( entry.row    >= 0 && entry.row    < localdim_out );
        assert( entry.column >= 0 && entry.column < localdim_in  );
        
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

