#ifndef INCLUDEGUARD_FEM_HEIGHTPRODUCTMATRIX
#define INCLUDEGUARD_FEM_HEIGHTPRODUCTMATRIX


#include <vector>
#include <iostream>
// #include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/densematrix.hpp"
#include "../operators/linearoperator.hpp"





inline DenseMatrix gradientproductmatrix( int n, int ambientdim, const std::vector<std::vector<Float>> list_of_vertices )
{
    
    assert( n <= ambientdim );
    assert( list_of_vertices.size() == n+1 )
    for( int v = 1; v <= n, v++ ) assert( list_of_vertices[v].size() == ambientdim );
    
    // create the matrix with the differences
    DenseMatrix Jacobian( ambientdim, n );
    
    for( int v = 1; v <= n, v++ )
    for( int c = 0; c < ambientdim; c++ )
        Jacobian.set( c, v, list_of_vertices[v][c] - list_of_vertices[0][c] );
    
    return gradientproductmatrix( n, ambientdim, Jacobian );
    
}




inline DenseMatrix gradientproductmatrix( int n, int ambientdim, DenseMatrix Jacobian )
{
    
    assert( n <= ambientdim );
    assert( Jacobian.getdimin() == n && Jacobian.getdimout() == ambientdim ); 
        
    DenseMatrix Product = Inverse( Transpose( Jacobian ) * Jacobian );
    
    DenseMatrix special( n, n+1, 0. );
    
    for( int i = 0; i < n; i++ ) { 
        special( i, i ) =  0.;
        special( i, n ) = -1.;
    }
    
    return special * Product * Transpose( special );
    
}



#endif
