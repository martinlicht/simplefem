#ifndef INCLUDEGUARD_FEM_FEECLAGRANGEINCLUSIONMATRIX
#define INCLUDEGUARD_FEM_FEECLAGRANGEINCLUSIONMATRIX


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







inline SparseMatrix FEECLagrangeInclusionMatrix( Mesh& mesh, int n, int r )
{
    
    // check whether the parameters are right 
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( r >= 1 );
    assert( r == 1 );
    
    const int num_simplices = mesh.count_simplices( n );
    const int num_vertices  = mesh.count_simplices( 0 );
    
    const int dim_out = num_simplices * (n+1);
    const int dim_in  = num_vertices;
    
    const int num_entries = num_simplices * (n+1);
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s  = 0; s  <  num_simplices; s++  )
    for( int vi = 0; vi <= n;             vi++ )
    {
        
        int index_of_entry = s * (n+1) + vi; 
            
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = s * (n+1) + vi;
        entry.column = mesh.get_subsimplex( n, 0, s, vi );
        entry.value  = 1.0;
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}




#endif
