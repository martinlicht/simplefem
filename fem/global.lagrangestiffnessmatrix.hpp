#ifndef INCLUDEGUARD_FEM_LAGRANGESTIFFNESSMATRIX
#define INCLUDEGUARD_FEM_LAGRANGESTIFFNESSMATRIX


#include <vector>
#include <iostream>
// #include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../mesh/mesh.hpp"





inline SparseMatrix LagrangeStiffnessMatrix( Mesh& mesh, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( r == 1 ); // only lowest order for the time being
    
    // Auxiliary calculations and preparations
    
    int n = mesh.getinnerdimension();
    
    const int num_volumes = mesh.count_simplices( n );
        
    const int num_vertices = mesh.count_simplices( 0 );
    
    const int dim_in  = num_vertices;
    const int dim_out = num_vertices;
    
    
    // Set up sparse matrix
    
    SparseMatrix ret( dim_out, dim_in, num_volumes * (n+1)*(n+1) );
    
    
    // go over the n simplices and their k subsimplices
    for( int t  = 0; t  <  num_volumes; t++  )
    for( int v1 = 0; v1 <= n;           v1++ )
    for( int v2 = 0; v2 <= n;           v2++ )
    {
        
        SparseMatrix::MatrixEntry entry;
        
        int index_of_entry = t * (n+1)*(n+1) + v1 * (n+1) + v2;
        
        entry.row    = mesh.get_subsimplex( n, 0, t, v1 );
        
        entry.column = mesh.get_subsimplex( n, 0, t, v2 );
        
        entry.value  = 0.;
        
        DenseMatrix Jac = mesh.getTransformationJacobian( n, t );
        
        DenseMatrix GradProds = mesh.getGradientProductMatrix( n, t );
        
        entry.value = GradProds( v1, v2 ) * Determinant( Jac ) / factorial( n );
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}



#endif
