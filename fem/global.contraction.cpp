
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/cholesky.hpp"
#include "../dense/functions.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/matrixtensorproduct.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/local.polynomialmassmatrix.hpp"

#include "../fem/global.massmatrix.hpp"



SparseMatrix FEECBrokenContractionMatrix( const Mesh& mesh, int n, int k, int r , int l, int s, FloatVector field )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( r >= 0 );
    assert( s >= 0 );
    assert( k >= 0 && k <= n );
    assert( l >= 0 && l <= n );
    assert( r >= s && k >= l );

    assert( k == l ); // restricted special case for now...
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim_in  = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    
    const int localdim_out = binomial_integer( n+(r+s), n ) * binomial_integer( n+1, (k-l) );
    
    const auto mis_input  = generateMultiIndices( IndexRange(0,n), r   );
    const auto mis_factor = generateMultiIndices( IndexRange(0,n), s   );
    const auto mis_output = generateMultiIndices( IndexRange(0,n), r+s );

    std::vector<std::array<int,3>> couplings;

    for( int m_i = 0; m_i < mis_input.size();  m_i++ )
    for( int m_f = 0; m_f < mis_factor.size(); m_f++ )
    {
        const auto& mi_input  = mis_input[m_i];
        const auto& mi_factor = mis_factor[m_f];
        
        const auto mi_prod = mi_input + mi_factor;
        
        for( int m_o = 0; m_o < mis_output.size(); m_o++ )
        {
            if( mi_prod == mis_output[m_o] ) {
                couplings.push_back( {m_i,m_f,m_o} );
            }                
        }
    }

    LOG << "Couplings: " << couplings.size() << nl;
    for( const auto& a : couplings ) LOG << a[0] << space << a[1] << space << a[2] << nl;


    
    const int dim_in      = num_simplices * localdim_in;
    const int dim_out     = num_simplices * localdim_out;
    const int num_entries = num_simplices * couplings.size() * binomial_integer( n+1, k );

    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        DenseMatrix GPM    = mesh.getGradientProductMatrix( n, s );
        DenseMatrix formMM = SubdeterminantMatrix( GPM, k );

        assert( not formMM.iszero() );
        assert( formMM.issquare() );
        assert( formMM.issymmetric() );

        for( int coupling_index = 0; coupling_index < couplings.size();        coupling_index++ ) 
        for( int f_i            = 0;            f_i < binomial_integer(n+1,k); f_i++            )
        {

            int m_i = couplings[coupling_index][0];
            int m_f = couplings[coupling_index][1];
            int m_o = couplings[coupling_index][2];

            SparseMatrix::MatrixEntry entry;
            entry.row    = s * mis_output.size() + m_o;
            entry.column = s * mis_input.size() * binomial_integer(n+1,k)  + m_i * binomial_integer(n+1,k) + f_i;
            
            entry.value = 0.;
            for( int g_i = 0; g_i < binomial_integer(n+1,k); g_i++ )
                entry.value += field[ s * mis_factor.size() * binomial_integer(n+1,k)  + m_f * binomial_integer(n+1,k) + f_i ] * formMM( f_i, g_i );

            int index_of_entry = s * couplings.size() * binomial_integer( n+1, k ) + coupling_index * binomial_integer( n+1, k ) + f_i;
            
            ret.setentry( index_of_entry,  entry );
        }        
        
    }
    
    ret.isfinite();

    return ret;
}





