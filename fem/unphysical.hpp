#ifndef INCLUDEGUARD_FEM_UNPHYSICAL_HPP
#define INCLUDEGUARD_FEM_UNPHYSICAL_HPP

#include <algorithm>
#include <vector>

#include "../basic.hpp"
#include "../utility/random.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/linearoperator.hpp"
#include "../operators/simpleoperators.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/matrixtensorproduct.hpp"
#include "../dense/functions.hpp"
#include "../mesh/mesh.hpp"
#include "../fem/utilities.hpp"


/************************
****
****  Unphysical operations
****  
************************/



inline FloatVector FEECCanonicalizeBroken( int n, int k, int r, const FloatVector& vec )
{
    Assert( 0 <= n );
    Assert( 0 <= r );
    Assert( 0 <= k and k <= n );
    vec.check();

    if( k == 0 ) return vec;

    const int poly_size = binomial_integer( n+r, n );
    const int form_size = binomial_integer( n+1, k );

    const int space_size = poly_size * form_size;
    
    const int cellcount = vec.getdimension() / space_size;
    
    Assert( cellcount * space_size == vec.getdimension(), cellcount, space_size, vec.getdimension(), n, k, r, poly_size, form_size );
    Assert( 0 == vec.getdimension() % space_size );
    
    DenseMatrix Aux1( n+1, n+1, 0. );
    for( int i = 1; i <= n; i++ ) {
        Aux1(i,i) = 1.;
        Aux1(i,0) = -1.;
    }
    
    DenseMatrix Aux2 = SubdeterminantMatrix( Aux1, k );

    assert( Aux2.issquare() and Aux2.getdimout() == form_size );
    
    FloatVector ret( vec.getdimension(), 0. );
    
    for( int c = 0; c < cellcount; c++ )
    for( int p = 0; p < poly_size; p++ )
    for( int i = 0; i < form_size; i++ )
    for( int j = 0; j < form_size; j++ )
        ret[ c * space_size + p * form_size + i ] += Aux2(i,j) * vec[ c * space_size + p * form_size + j ];
    
    return ret;
}


inline FloatVector FEECRandomizeBroken( int n, int k, int r, const FloatVector& vec )
{
    Assert( 0 <= n );
    Assert( 0 <= r );
    Assert( 0 <= k and k <= n );
    vec.check();

    if( k == 0 ) return vec;

    const int poly_size = binomial_integer( n+r, n );
    const int form_size = binomial_integer( n+1, k );

    const int space_size = poly_size * form_size;
    
    const int cellcount = vec.getdimension() / space_size;
    
    Assert( cellcount * space_size == vec.getdimension(), cellcount, space_size, vec.getdimension(), n, k, r, poly_size, form_size );
    Assert( 0 == vec.getdimension() % space_size );
    
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
    
    DenseMatrix Aux2 = SubdeterminantMatrix( Aux1, k );

    assert( Aux2.issquare() and Aux2.getdimout() == form_size );
    
    FloatVector ret( vec.getdimension(), 0. );
    
    for( int c = 0; c < cellcount; c++ )
    for( int p = 0; p < poly_size; p++ )
    for( int i = 0; i < form_size; i++ )
    for( int j = 0; j < form_size; j++ )
        ret[ c * space_size + p * form_size + i ] += Aux2(i,j) * vec[ c * space_size + p * form_size + j ];
    
    return ret;
}
 

#endif
