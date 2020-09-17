#ifndef INCLUDEGUARD_FEM_LAGRANGEMASSMATRIX
#define INCLUDEGUARD_FEM_LAGRANGEMASSMATRIX


// #include <cassert>
#include <iostream>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../mesh/mesh.hpp"




///////////////////////////////////////////////////////
//                                                   //
//  Matrix for the continuous Lagrange mass pairing  //
//  with poly degree r                               //
//                                                   //
///////////////////////////////////////////////////////


inline SparseMatrix LagrangeMassMatrix( const Mesh& mesh, int r )
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
        
        int vertex1 = mesh.get_subsimplex( n, 0, t, v1 );
        int vertex2 = mesh.get_subsimplex( n, 0, t, v2 );
        
        entry.row    = mesh.get_subsimplex( n, 0, t, v1 );
        entry.column = mesh.get_subsimplex( n, 0, t, v2 );
        
//         DenseMatrix Jac = mesh.getTransformationJacobian( n, t );
        Float measure = mesh.getMeasure( n, t );
        
        if( v1 == v2 )
            entry.value = 2. * factorial_numerical(n) * measure / factorial_numerical( 2 + n );
        else
            entry.value =      factorial_numerical(n) * measure / factorial_numerical( 2 + n );
                
        if( mesh.get_flag( 0, vertex1 ) == SimplexFlagDirichlet or mesh.get_flag( 0, vertex2 ) == SimplexFlagDirichlet )
            entry.value = 0.;

        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}



#endif
