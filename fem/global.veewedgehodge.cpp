
#include <cmath>

#include <algorithm>
#include <array>
#include <vector>

#include "../base/include.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"
#include "../fem/polynomialmassmatrix.hpp"

#include "../fem/global.veewedgehodge.hpp"





// =============================================================================
// Generates a vector whose product with any scalar finite element field
// produces the integral of that field.
// =============================================================================
FloatVector FEECScalarIntegral( const Mesh& mesh, int n, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    
    // Auxiliary calculations and preparations
    
	const int num_simplices = mesh.count_simplices( n );
    
    const int localdim = binomial_integer( n + r, n );
    
    const int num_entries = num_simplices * localdim;
    
    FloatVector ret( num_entries );
	
	// Numerical data : for the standard volume form
    
    FloatVector polynomial_weights( binomial_integer( n + r, n ) );
	
	std::vector<MultiIndex> multis = generateMultiIndices( IndexRange(0,n), r );
    
    const int N = multis.size();
    
    assert( N == binomial_integer( n + r, r ) );
    
    for( int index_poly = 0; index_poly < N; index_poly++ )
    {
        MultiIndex alpha = multis[index_poly];
        assert( absolute( alpha ) == r ); 
        // n! alpha! / (n+|alpha|)!
        polynomial_weights[index_poly] = factorial_numerical(n) * alpha.factorial_numerical() / factorial_numerical( n + r ); 
    }
	
	if( r == 0 )
    for( auto w : polynomial_weights ) Assert( is_numerically_one( w ) );
	
	// Fill in the data into the return vector 
	for( int s = 0;          s < num_simplices; s++          )
	for( int index_poly = 0; index_poly < N;    index_poly++ )
	{
		Float volume = mesh.getMeasure( n, s );
        ret[ s * N + index_poly] = volume * polynomial_weights[index_poly];
	}

    for( auto v : ret ) assert( v > 0. );
	
	return ret; 
}






// =============================================================================
// Generates a vector whose product with any volume finite element field
// produces the integral of that field.
// =============================================================================
FloatVector FEECVolumeFormIntegral( const Mesh& mesh, int n, int r )
{
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
    const int localdim      = binomial_integer( n + r, n ) * (n+1);
    const int num_entries   = num_simplices * localdim;
    
    FloatVector ret( num_entries );
	
	// Numerical data : for the standard volume form
    
    std::vector<MultiIndex> multis = generateMultiIndices( IndexRange(0,n), r );
    
    const int N = multis.size();
    
    assert( N == binomial_integer( n + r, r ) );
    
    FloatVector polynomial_weights( N );
	
	for( int index_poly = 0; index_poly < N; index_poly++ )
    {
        MultiIndex alpha = multis[index_poly];
        assert( absolute( alpha ) == r ); // n! alpha! / (n+|alpha|)!
        polynomial_weights[index_poly] = factorial_numerical(n) * alpha.factorial_numerical() / factorial_numerical( n + r ); 
    }

    // Numerical data: catalog the signs
	
	std::vector<IndexMap> sigmas = generateSigmas( IndexRange(1,n), IndexRange(0,n) );
    assert( sigmas.size() == n+1 );
	
	FloatVector signs( sigmas.size() );
	
    for( int index_form = 0; index_form < sigmas.size(); index_form++ )
	for( int p = 0;          p <= n;                     p++          )
	{
		const auto& sigma = sigmas[index_form];
		
		bool p_found = sigma.has_value_in_range(p);

        if( p_found ) continue;
		
        assert( not std::isfinite( signs[index_form] ) );

        signs[index_form] = sign_power( p );

        if( p == 0 ) assert( sign_power(p) ==  1 );
        if( p == 1 ) assert( sign_power(p) == -1 );
        if( p == 2 ) assert( sign_power(p) ==  1 );
	}
    
	// Fill in the data into the return vector 
	
    for( int s = 0;          s < num_simplices;          s++          )
    for( int index_poly = 0; index_poly < N;             index_poly++ )
	for( int index_form = 0; index_form < sigmas.size(); index_form++ )
	{
		auto Jac = mesh.getTransformationJacobian(n,s);
        assert( Jac.is_square() );
        auto orientation = sign_integer( Determinant(Jac) );

        // LOG << s << space << orientation << space << signs[index_form] << nl;
        
        Float value = orientation * signs[index_form] * polynomial_weights[index_poly] / factorial_numerical( n );

        ret[ s * N * (n+1) + index_poly * (n+1) + index_form] = value; 
	}


	
	return ret; 
}













// =============================================================================
// Vee product of fixed field of poly_s l-forms with poly_r k-forms 
// 
// Currently only pairings that lead to scalar fields
// =============================================================================
SparseMatrix FEECBrokenVeeMatrix( const Mesh& mesh, int n, int form_k, int poly_r , int form_l, int poly_s, const FloatVector& field )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( poly_r >= 0 );
    assert( poly_s >= 0 );
    assert( form_k >= 0 && form_k <= n );
    assert( form_l >= 0 && form_l <= n );
    assert( form_k >= form_l );
    assert( field.is_finite() );

    assert( form_k == form_l ); // restricted special case for now. Requires understanding of pairing between general forms
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim_in  = binomial_integer( n + ( poly_r          ), n ) * binomial_integer( n+1, ( form_k        ) );
    const int localdim_out = binomial_integer( n + ( poly_r + poly_s ), n ) * binomial_integer( n+1, ( form_k-form_l ) );
    
    const auto mis_input  = generateMultiIndices( IndexRange(0,n), poly_r          );
    const auto mis_factor = generateMultiIndices( IndexRange(0,n),          poly_s );
    const auto mis_output = generateMultiIndices( IndexRange(0,n), poly_r + poly_s );

    // We save the products of multiindices 

    std::vector<std::array<int,3>> mi_couplings;

    for( int m_i = 0; m_i < mis_input.size();  m_i++ )
    for( int m_f = 0; m_f < mis_factor.size(); m_f++ )
    {
        const auto& mi_input  = mis_input[m_i];
        const auto& mi_factor = mis_factor[m_f];
        
        const auto mi_prod = mi_input + mi_factor;
        
        int m_o = 0;

        for(; m_o < mis_output.size() && mi_prod != mis_output[m_o]; m_o++ );

        assert( m_o != mis_output.size() );
        assert( mi_prod == mis_output[m_o] );

        mi_couplings.push_back( {m_i,m_f,m_o} );
    }

    const int dim_in      = num_simplices * localdim_in;
    const int dim_out     = num_simplices * localdim_out;
    const int num_entries = num_simplices * mi_couplings.size() * binomial_integer( n+1, form_k );

    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        // we compute the local product matrix for the form part

        DenseMatrix GPM    = mesh.getGradientProductMatrix( n, s );
        DenseMatrix formMM = SubdeterminantMatrix( GPM, form_k );
        
        assert( formMM.is_finite() );
        assert( not formMM.is_zero() );
        assert( formMM.is_square() );
        assert( formMM.is_symmetric() );

        // Iterate over multiindex couplings and input form degrees.
        for( int coupling_index = 0; coupling_index < mi_couplings.size();          coupling_index++ ) 
        for( int f_i            = 0;            f_i < binomial_integer(n+1,form_k); f_i++            )
        {

            int m_i = mi_couplings[coupling_index][0]; // input 
            int m_f = mi_couplings[coupling_index][1]; // output 
            int m_o = mi_couplings[coupling_index][2]; // pairing

            int row    = s * mis_output.size() * 1                             + m_o;
            int column = s * mis_input.size()  * binomial_integer(n+1,form_k)  + m_i * binomial_integer(n+1,form_k) + f_i;
            
            Float value = 0.;
            for( int g_i = 0; g_i < binomial_integer(n+1,form_l); g_i++ )
                value += field[ s * mis_factor.size() * binomial_integer(n+1,form_l)  + m_f * binomial_integer(n+1,form_l) + g_i ] * formMM( g_i, f_i );

            Assert( std::isfinite(value) );

            int index_of_entry = s * mi_couplings.size() * binomial_integer( n+1, form_k ) + coupling_index * binomial_integer( n+1, form_k ) + f_i;
            
            ret.setentry( index_of_entry, { row, column, value } );
        }        
        
    }
    
    assert( ret.is_finite() );

    return ret;
}












// =============================================================================
// Wedge product of fixed field of poly_s l-forms with poly_r k-forms 
// =============================================================================

SparseMatrix FEECBrokenWedgeMatrix( const Mesh& mesh, int n, int form_k, int poly_r , int form_l, int poly_s, const FloatVector& field )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( poly_r >= 0 );
    assert( poly_s >= 0 );
    assert( form_k >= 0 && form_k <= n );
    assert( form_l >= 0 && form_l <= n );
    assert( field.is_finite() );

    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim_in  = binomial_integer( n + ( poly_r          ), n ) * binomial_integer( n+1, form_k          );    
    const int localdim_out = binomial_integer( n + ( poly_r + poly_s ), n ) * binomial_integer( n+1, form_k + form_l );
    
    // We save the products of multiindices 

    const auto mis_input  = generateMultiIndices( IndexRange(0,n), poly_r          );
    const auto mis_factor = generateMultiIndices( IndexRange(0,n), poly_s          );
    const auto mis_output = generateMultiIndices( IndexRange(0,n), poly_r + poly_s );

    std::vector<std::array<int,3>> mi_couplings;

    for( int m_i = 0; m_i < mis_input.size();  m_i++ )
    for( int m_f = 0; m_f < mis_factor.size(); m_f++ )
    {
        const auto& mi_input  = mis_input[m_i];
        const auto& mi_factor = mis_factor[m_f];
        
        const auto mi_prod = mi_input + mi_factor;
        
        int m_o = 0;

        for(; m_o < mis_output.size() && mi_prod != mis_output[m_o]; m_o++ );

        assert( m_o != mis_output.size() );
        assert( mi_prod == mis_output[m_o] );

        mi_couplings.push_back( {m_i,m_f,m_o} );
    }


    // We save the wedge of l-forms (factor) with k-forms (input) to get k+l-forms (output)

    const auto sigmas_input  = generateSigmas( IndexRange(1, form_k          ), IndexRange(0,n) );
    const auto sigmas_factor = generateSigmas( IndexRange(1,          form_l ), IndexRange(0,n) );
    const auto sigmas_output = generateSigmas( IndexRange(1, form_k + form_l ), IndexRange(0,n) );

    std::vector<std::array<int,4>> sigma_couplings;

    for( int s_i = 0; s_i < sigmas_input.size();  s_i++ )
    for( int s_f = 0; s_f < sigmas_factor.size(); s_f++ )
    {
        const auto& sigma_input  = sigmas_input[s_i];
        const auto& sigma_factor = sigmas_factor[s_f];
        
        int signum;
        IndexMap sigma_prod = mergeSigmas( sigma_factor, sigma_input, signum );

        if( signum == 0 ) continue;
        
        // for( int s_o = 0; s_o < sigmas_output.size(); s_o++ )
        // {
        //     if( sigma_prod == sigmas_output[s_o] ) {
        //         sigma_couplings.push_back( {s_i,s_f,s_o,signum} );
        //     }                
        // }

        int s_o = 0;

        for(; s_o < sigmas_output.size() && sigma_prod != sigmas_output[s_o]; s_o++ );

        assert( s_o != sigmas_output.size() );
        assert( sigma_prod == sigmas_output[s_o] );

        sigma_couplings.push_back( { s_i, s_f, s_o, signum } );
    }

    
    
    // prepare the final sparse matrix 

    const int dim_in      = num_simplices * localdim_in;
    const int dim_out     = num_simplices * localdim_out;
    const int num_entries = num_simplices * mi_couplings.size() * sigma_couplings.size();

    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        // For each simplex, we run over all possible product terms 
        // This variable saves how many terms have been processed:
        int local_index_of_entry = 0;

        // All multiindex and sigma couplings ...
        for( int mi_coupling_index = 0;    mi_coupling_index    < mi_couplings.size();    mi_coupling_index++    ) 
        for( int sigma_coupling_index = 0; sigma_coupling_index < sigma_couplings.size(); sigma_coupling_index++ ) 
        {

            int m_i = mi_couplings[mi_coupling_index][0];
            int m_f = mi_couplings[mi_coupling_index][1];
            int m_o = mi_couplings[mi_coupling_index][2];

            int s_i = sigma_couplings[sigma_coupling_index][0];
            int s_f = sigma_couplings[sigma_coupling_index][1];
            int s_o = sigma_couplings[sigma_coupling_index][2];

            int row     =                 s * mis_output.size() * binomial_integer( n+1, form_k + form_l ) + m_o * binomial_integer( n+1, form_k + form_l ) + s_o;
            int column  =                 s * mis_input.size()  * binomial_integer( n+1, form_k          ) + m_i * binomial_integer( n+1, form_k          ) + s_i;
            
            auto signum      = sigma_couplings[sigma_coupling_index][3];
            auto orientation = mesh.getOrientation(s);

            Float value = signum * orientation * field[ s * mis_factor.size() * binomial_integer( n+1,          form_l ) + m_f * binomial_integer( n+1,          form_l ) + s_f ];
            
            assert( std::isfinite(value) );

            int index_of_entry = s * mi_couplings.size() * sigma_couplings.size() + local_index_of_entry;
            
            ret.setentry( index_of_entry, { row, column, value } );
            local_index_of_entry++;
        }        
        
        assert( local_index_of_entry == mi_couplings.size() * sigma_couplings.size() );
    }
    
    assert( ret.is_finite() );

    return ret;
}













// =============================================================================
// Hodge Star Matrix (exported)
// =============================================================================

SparseMatrix FEECBrokenHodgeStarMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    
    /* 0. Prepare the auxiliary data */

    const int num_simplices = mesh.count_simplices( n );

    // input k forms 
    auto sigmas_input  = generateSigmas( IndexRange(1,k  ), IndexRange(0,n) );
    
    // test k forms, we only use those which do not start at 0
    auto sigmas_test   = generateSigmas( IndexRange(1,k  ), IndexRange(0,n) );
    {
        auto sigma_test_startnonzero = std::remove_if( sigmas_test.begin(), sigmas_test.end(),
                        [k](IndexMap im) -> bool { assert( im.is_strictly_ascending() ); return   k != 0 and im[1] == 0; }
        );
        sigmas_test.erase( sigma_test_startnonzero, sigmas_test.end() );
    }

    // output n-k forms 
    auto sigmas_output = generateSigmas( IndexRange(1,n-k), IndexRange(0,n) );
    
    // actual n-k forms for output, we only use those which do not start at 0
    auto sigmas_hodge  = generateSigmas( IndexRange(1,n-k), IndexRange(0,n) );
    {
        auto sigma_hodge_startnonzero = std::remove_if( sigmas_hodge.begin(), sigmas_hodge.end(),
                [n,k](IndexMap im) -> bool { assert( im.is_strictly_ascending() ); return n-k != 0 and im[1] == 0; }
        );
        sigmas_hodge.erase( sigma_hodge_startnonzero, sigmas_hodge.end() );
    }

    const int dim_polynomials = binomial_integer( n+r, n   );
    const int dim_sigmas_in   = binomial_integer( n+1,   k );
    const int dim_sigmas_out  = binomial_integer( n+1, n-k );
    
    assert( sigmas_input.size()  == dim_sigmas_in       );
    assert( sigmas_output.size() == dim_sigmas_out      );
    Assert( sigmas_test.size()   == sigmas_hodge.size(), n, k, sigmas_test.size(), sigmas_hodge.size() );
    
    
    /* 1. Assemble the algebraic matrix over the test and hodge forms. */ 
    /*    This corresponds to an Euclidean Hodge star.                 */ 

    const IndexMap standard_volume_form( IndexRange(1,n), IndexRange(0,n), [](int i)->int{ return i; } );

    DenseMatrix wedge_matrix( sigmas_test.size(), sigmas_hodge.size(), 0. );

    for( int s_h = 0; s_h < sigmas_hodge.size(); s_h++ )
    for( int s_t = 0; s_t < sigmas_test.size();  s_t++ )
    {
        const auto& sigma_hodge = sigmas_hodge[s_h];
        const auto& sigma_test  = sigmas_test[s_t];
        
        int signum;
        IndexMap sigma_prod = mergeSigmas( sigma_test, sigma_hodge, signum );

        if( signum == 0 ) continue;

        assert( sigma_prod.is_strictly_ascending() and sigma_prod == standard_volume_form );

        wedge_matrix( s_t, s_h ) = signum;
    }

    assert( wedge_matrix.is_square() );

    const DenseMatrix wedge_matrix_inv = Inverse( wedge_matrix ); // TODO(martinlicht): This inversion can be made much simpler ... 

    assert( ( wedge_matrix_inv * wedge_matrix - IdentityMatrix( wedge_matrix.getdimin() ) ).is_numerically_small() );
    assert( ( Transpose(wedge_matrix) * wedge_matrix - IdentityMatrix( wedge_matrix.getdimin() ) ).is_numerically_small() );

    
    /* 2. Extend this matrix to the full set of indices, including those that start with zero. */ 

    DenseMatrix wedge_matrix_inv_full( sigmas_output.size(), sigmas_input.size(), 0. );
    
    for( int s_h = 0; s_h < sigmas_hodge.size();  s_h++ )
    for( int s_t = 0; s_t < sigmas_test.size();   s_t++ )
    for( int s_i = 0; s_i < sigmas_input.size();  s_i++ )
    for( int s_o = 0; s_o < sigmas_output.size(); s_o++ )
    {
        if( sigmas_test[s_t] != sigmas_input[s_i] or sigmas_hodge[s_h] != sigmas_output[s_o] ) {
            continue;
        } else {
            wedge_matrix_inv_full( s_o, s_i ) = wedge_matrix_inv( s_h, s_t );
        }
    }

    assert( wedge_matrix_inv_full.is_finite() );
    // for( auto x : sigmas_output ) LOG << x << nl;
    // for( auto x : sigmas_input ) LOG << x << nl;
    // LOG << nl << wedge_matrix_inv_full << nl;
    
    
    
    
    /* 3. Run over the simplices and compile the local entries. */
    
    /*                                                          */
    /* The Hodge star is locally characterized by the equation  */
    /* intvolume( v wedge star(u) ) = int mass( v, u )          */
    /* where                                                    */
    /*                                                          */
    /* - u is the input form (input sigmas)                     */
    /* - star(u) is the Hodge form (hodge sigmas)               */
    /* - v is some test form (test forms)                       */
    /*                                                          */
    /* In terms of matrices, H star(u) = M u, which means       */
    /* star(u) = inv(H) M u                                     */
    /*                                                          */
    /* The last inverse might be understood as a pseudoinverse. */
    /*                                                          */
    
    const int dim_in      = num_simplices * dim_polynomials * dim_sigmas_in;
    const int dim_out     = num_simplices * dim_polynomials * dim_sigmas_out;
    const int num_entries = num_simplices * dim_polynomials * dim_sigmas_in * dim_sigmas_out;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        /* 3.1 Build the local mass matrix. */
        
        const DenseMatrix GPM               = mesh.getGradientProductMatrix( n, s );
        
        const DenseMatrix formMM            = SubdeterminantMatrix( GPM, k ) * mesh.getMeasure( n, s ) * factorial_integer(n);

        assert( formMM.is_finite() );
        
        const DenseMatrix full_local_matrix = wedge_matrix_inv_full * formMM;

        /* 3.1 Get the full matrix. */
        
        assert( full_local_matrix.is_finite() );
        
        const Float scaling = 1.; 
    
        for( int p   = 0; p   < dim_polynomials; p++   )
        for( int s_i = 0; s_i < dim_sigmas_in;   s_i++ )
        for( int s_o = 0; s_o < dim_sigmas_out;  s_o++ )
        {

            int index_of_entry = s * dim_polynomials * dim_sigmas_in * dim_sigmas_out + p * dim_sigmas_in * dim_sigmas_out + s_i * dim_sigmas_out + s_o;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * dim_polynomials * dim_sigmas_out + p * dim_sigmas_out + s_o;
            entry.column = s * dim_polynomials * dim_sigmas_in  + p * dim_sigmas_in  + s_i;
            entry.value  = full_local_matrix( s_o, s_i ) * scaling;
                
            ret.setentry( index_of_entry, entry );

        }
        
    }
    
    return ret;
}





















// =============================================================================
// Provides the Euclidean Hodge Star Matrix
// =============================================================================

DenseMatrix EuclideanHodgeStar( int n, int k )
{
    assert( n >= 0 );
    assert( k >= 0 && k <= n );
    
    const std::vector<IndexMap> sigmas_in  = generateSigmas( IndexRange( 1, k   ), IndexRange( 1, n ) );
    const std::vector<IndexMap> sigmas_out = generateSigmas( IndexRange( 1, n-k ), IndexRange( 1, n ) );

    assert( sigmas_in.size() == sigmas_out.size() );

    DenseMatrix ret( sigmas_in.size(), sigmas_out.size(), 0. );

    for( int r = 0; r < sigmas_out.size(); r++ )
    {
        const auto& sigma_out = sigmas_out[r];

        int c = 0;

        for( ; c < sigmas_in.size(); c++ ) 
        {
            const auto& sigma_in = sigmas_in[c];

            bool common_index = false;

            for( auto v : sigma_out.getvalues() )
            for( auto w :  sigma_in.getvalues() )
            {
                common_index = common_index or ( v == w );
            }

            if( not common_index ) break;
        }

        assert( c != sigmas_in.size() );

        const auto& sigma_in = sigmas_in[c];

        int signum = 1.; // determine the specific coefficient needed

        for( int i = 1; i <= k;   i++ )
        for( int j = 1; j <= n-k; j++ )
        {
            if( sigma_in[i] > sigma_out[j] ) signum *= -1;
        }

        ret(r,c) = signum;
    }

    return ret;
}




















// =============================================================================
// Hodge Star Matrix (alternative implementation)
// =============================================================================

SparseMatrix FEECBrokenHodgeStarMatrix_Alternative( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    
    const int dim_in      = num_simplices * localdim;
    const int dim_out     = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim;
    
    const std::vector<MultiIndex> multis_in  = generateMultiIndices( IndexRange( 0, n ), r );
    const std::vector<IndexMap>   sigmas_in  = generateSigmas( IndexRange( 1, k   ), IndexRange( 0, n ) );
    
    const std::vector<MultiIndex> multis_out = generateMultiIndices( IndexRange( 0, n ), r );
    const std::vector<IndexMap>   sigmas_out = generateSigmas( IndexRange( 1, n-k ), IndexRange( 0, n ) );
    
    // TODO(martinlicht): complete code 
    // assert( multis_dest.size() * sigmas_dest.size() == multis_src.size()  * sigmas_src.size() );
    assert( dim_out == dim_in );
    
    
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    const DenseMatrix polyMM = polynomialmassmatrix( n, r );

    assert( polyMM.is_square() and polyMM.getdimin() == binomial_integer( n+r, n ) );
    
    /*
    Method 1:

    1. Develop a QR decomposition of the transformation matrix to obtain local Euclidean coordinate system;
       In the full-dimensional case, an alternative method is feasible as well (introduce boolean switch).
    2. Write barycentric derivatives in terms of those coordinates (canonicalize the zeroth gradient)
    3. Apply Euclidean Hodge star
    4. Transform back into barycentric coordinate system 

    Method 2:
    1. Write the mass pairing of the n-k forms 
    2. Define the Hodge star implicitly.
    */

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        // obtain geometric information 
        
        Float measure      = mesh.getMeasure( n, s );

        DenseMatrix GPM    = mesh.getGradientProductMatrix( n, s );
          
        assert( measure >= 0. );

        // create local mass matrix 
        
        DenseMatrix formMM = SubdeterminantMatrix( GPM, k );
    
        DenseMatrix fullMM = MatrixTensorProduct( polyMM, formMM ) * measure;

        // create local Hodge product matrix 
        
        DenseMatrix formHM(0); // TODO(martinlicht): wierd ...

        DenseMatrix fullHM = MatrixTensorProduct( polyMM, formHM ) * measure;

        // compute the final local matrix and register its entries 

        DenseMatrix localmatrix = Inverse(fullMM) * fullHM;

        for( int i = 0; i < localdim; i++ )
        for( int j = 0; j < localdim; j++ )
        {
            int index_of_entry = s * localdim * localdim + i * localdim + j;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * localdim + i;
            entry.column = s * localdim + j;
            entry.value  = localmatrix( i, j );
            
            ret.setentry( index_of_entry, entry );
        }
        
    }
    
    return ret;
}




