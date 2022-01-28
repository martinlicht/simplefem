
#include <iostream>
#include <vector>

#include "../basic.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/lagrangematrices.hpp"



SparseMatrix LagrangeBrokenMassMatrix( const Mesh& mesh, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( r == 1 ); // only lowest order for the time being
    
    // Auxiliary calculations and preparations
    
    int n = mesh.getinnerdimension();
    
    const int num_volumes = mesh.count_simplices( n );
    
    const int dim_in  = num_volumes * (n+1);
    const int dim_out = num_volumes * (n+1);
    
    
    // Set up sparse matrix
    
    SparseMatrix ret( dim_out, dim_in, num_volumes * (n+1)*(n+1) );
    
    
    // go over the n simplices and their k subsimplices
    for( int t  = 0; t  <  num_volumes; t++  )
    for( int v1 = 0; v1 <= n;           v1++ )
    for( int v2 = 0; v2 <= n;           v2++ )
    {
        
        SparseMatrix::MatrixEntry entry;
        
        int index_of_entry = t * (n+1)*(n+1) + v1 * (n+1) + v2;
        
        entry.row    = t * (n+1) + v1;
        
        entry.column = t * (n+1) + v2;
        
        Float measure = mesh.getMeasure( n, t );
        
        if( v1 == v2 )
            entry.value = 2. * factorial_numerical(n) * measure / factorial_numerical( 2 + n );
        else
            entry.value =      factorial_numerical(n) * measure / factorial_numerical( 2 + n );
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}




SparseMatrix LagrangeMassMatrix( const Mesh& mesh, int r )
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




SparseMatrix LagrangeBrokenStiffnessMatrix( const Mesh& mesh, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( r == 1 ); // only lowest order for the time being
    
    // Auxiliary calculations and preparations
    
    int n = mesh.getinnerdimension();
    
    const int num_volumes = mesh.count_simplices( n );
    
    const int dim_in  = num_volumes * (n+1);
    const int dim_out = num_volumes * (n+1);
    
    
    // Set up sparse matrix
    
    SparseMatrix ret( dim_out, dim_in, num_volumes * (n+1)*(n+1) );
    
    
    // go over the n simplices and their k subsimplices
    for( int t  = 0; t  <  num_volumes; t++  )
    for( int v1 = 0; v1 <= n;           v1++ )
    for( int v2 = 0; v2 <= n;           v2++ )
    {
        
        SparseMatrix::MatrixEntry entry;
        
        int index_of_entry = t * (n+1)*(n+1) + v1 * (n+1) + v2;
        
        entry.row    = t * (n+1) + v1;
        
        entry.column = t * (n+1) + v2;
        
//         entry.value  = 0.;
        
        DenseMatrix Jac = mesh.getTransformationJacobian( n, t );
        
        DenseMatrix GradProds = mesh.getGradientProductMatrix( n, t );
        
        Float measure = mesh.getMeasure( n, t );
        
        entry.value = GradProds( v1, v2 ) * measure; // /factorial_numerical( n );
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}





SparseMatrix LagrangeStiffnessMatrix( const Mesh& mesh, int r )
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

        entry.row    = vertex1; 
        entry.column = vertex2;
        
        //DenseMatrix Jac = mesh.getTransformationJacobian( n, t );
        
        DenseMatrix GradProds = mesh.getGradientProductMatrix( n, t );
        
        Float measure = mesh.getMeasure( n, t );
        
        entry.value = GradProds( v1, v2 ) * measure; // / factorial_numerical( n );
        
        if( mesh.get_flag( 0, vertex1 ) == SimplexFlagDirichlet or mesh.get_flag( 0, vertex2 ) == SimplexFlagDirichlet )
            entry.value = 0.;
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}





SparseMatrix LagrangeInclusionMatrix( const Mesh& mesh, int n, int r )
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
        
        int vertex = mesh.get_subsimplex( n, 0, s, vi );
        
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = s * (n+1) + vi;
        entry.column = vertex;
        entry.value  = 1.0;
        
        if( mesh.get_flag( 0, vertex ) == SimplexFlagDirichlet )
            entry.value = 0.0;
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}



