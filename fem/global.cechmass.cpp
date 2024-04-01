
#include <vector>

#include "../basic.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/global.cechmass.hpp"

SparseMatrix FEECBrokenCechMatrix( const Mesh& mesh, int n, int k, int s )
{
    
    assert( s >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
    
    const int num_chains    = mesh.count_simplices( k );
    
    const int chains_per_simplex = binomial_integer( n + 1, k + 1 );
    
    const int dim_in      = num_chains;
    const int dim_out     = num_chains;
    const int num_entries = num_simplices * chains_per_simplex;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices;      s++ )
    for( int i = 0; i < chains_per_simplex; i++ )
    {
        
        const Float measure = mesh.getMeasure( n, s );
        assert( measure >= 0. );

        const int index = mesh.get_subsimplex( n, k, s, i );

        Float diameter = 0.;

        if( k > 0 ) {
            diameter = mesh.getDiameter( k, index );
        } else {
            auto edges = mesh.getsupersimplices( 1, 0, index );
            for( int e : edges )
                diameter = maximum( diameter, mesh.getDiameter(1,e) );
        }
        assert( diameter > 0 );

        int index_of_entry = s * chains_per_simplex + i;
            
        SparseMatrix::MatrixEntry entry;
        entry.row    = index;
        entry.column = index;
        entry.value  = measure * power_numerical( diameter, n - k );
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}



