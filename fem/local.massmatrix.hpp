#ifndef INCLUDEGUARD_FEM_MASSMATRIX
#define INCLUDEGUARD_FEM_MASSMATRIX


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





inline DenseMatrix elementmassmatrix( int n, int ambientdim, int r, int k, DenseMatrix Jacobian )
{
    
    assert( Jacobian.getdimin() == n && Jacobian.getdimout() == ambientdim ); 
    assert( n <= ambientdim );
    assert( n >= 0 && n >= r );
        
    DenseMatrix polyMM = polynomialmassmatrix( n, r );
    
    DenseMatrix gradMM = gradientproductmatrix( n, ambientdim, Jacobian );
    
    DenseMatrix formMM = SubdeterminantMatrix( gradMM );
    
    return TensorProduct( polyMM, formMM ) * determinant( Jacobian ) / factorial( n );
    
}



#endif
