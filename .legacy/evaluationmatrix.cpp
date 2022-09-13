// 
// Suppose that the columns of lpsbc are the interpolation points
// (in barycentric coordinates) to interpolate degree r polynomials 
// over a d simplex exactly.
// 
// Suppose that mis cointains a collection of multiindices
// over a d simplex 
// 
// The output matrix is as follows:
// - columns correspond to the interpolation points 
// - rows correspond to the multiindices (Lagrange basis)
// - the entries are value of the corresponding polynomial at the corresponding point
// 
// TODO Reread and compare with previous version 

inline DenseMatrix PointValuesOfMonomials( std::vector<MultiIndex> mis, const DenseMatrix& lpsbc )
{
    
    const int num_points = lpsbc.getdimin();
    const int num_mis    = mis.size();    
    const int dim        = lpsbc.getdimout() - 1;
        
    DenseMatrix ret( num_points, num_mis );
    
    for( int c = 0; c < num_mis;    c++ ) // c -> barycentric poly 
    for( int r = 0; r < num_points; r++ ) // r -> interpolation point
    {
        assert( mis[c].getSourceRange().max() == dim+1 );
        assert( mis[c].getSourceRange().min() == 0     );

        ret(r,c) = 1.;
        for( int d = 0; d <= dim; d++ )
            if( mis[c][d] != 0 )
                ret(r,c) *= power_numerical( lpsbc(d,r), mis[c][d] );
    }
    
    return ret;
    
}






















// Given the Jacobian of the coordinate transformation, 
// this produces a matrix that transform from Euclidean coordinates 
// to barycentric coordinates 

// inline DenseMatrix BarycentricProjectionMatrixALTERVERSUCH( const DenseMatrix& J )
// {
//     assert( J.getdimout() >= J.getdimin() );
    
//     // We implement J^+ = inv( J^t J ) J^t
    
//     const auto Jt = Transpose(J);
    
//     const auto F = Inverse( Jt * J ) * Jt; // n x d
    
//     DenseMatrix ret( J.getdimin()+1, J.getdimout(), 0.0 );
    
//     for( int r = 0; r < J.getdimin();  r++ )
//     for( int c = 0; c < J.getdimout(); c++ )
//         ret( r+1, c ) = F(r,c);
    
//     return ret;
// }















