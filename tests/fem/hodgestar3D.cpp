

/**/

#include <cmath>

#include <string>

#include "../../base/include.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.hodgestarmatrix.hpp"
#include "../../utility/convergencetable.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: (3D) Hodge star behaves as intended" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = UnitCube3D();

    M.check();
    
    

    /* test the Euclidean Hodge matrix */
    for( int n = 1; n <= 5; n++ )
    for( int k = 0; k <= n; k++ )
    {
        auto hodge = EuclideanHodgeStar( n,   k );
        auto other = EuclideanHodgeStar( n, n-k );

        auto prod = other * hodge;

        auto signed_identity = sign_power( k * (n-k) ) * IdentityMatrix( other.getdimout() );

        LOG << prod << nl << signed_identity << nl;
        Assert( ( prod - signed_identity ).norm() == 0., prod, signed_identity );
    }
    

    
    
    const int r_min = 0;
    
    const int r_max = 1;
    
    const int l_min = 0;
    
    const int l_max = 3;

    const int number_of_random_trial_samples = 10;
    
    const int number_of_random_test_samples  = 10;
    
    Float inv_errors_scalar[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float inv_errors_vector[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float inv_errors_pseudo[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float inv_errors_volume[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    
    Float mass_errors_scalar[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float mass_errors_vector[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float mass_errors_pseudo[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float mass_errors_volume[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    
    
    
    for( int l = 0; l < l_min; l++ )
        M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ ){
        
        LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
        
        for( int r = r_min; r <= r_max; r++ ) 
        {
            
            LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;

            LOG << "assemble mass matrices..." << nl;
            
            SparseMatrix massmatrix_scalar = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
            
            SparseMatrix massmatrix_vector = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
            
            SparseMatrix massmatrix_pseudo = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );
            
            SparseMatrix massmatrix_volume = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r );
            
            assert( massmatrix_scalar.is_finite() );
            assert( massmatrix_vector.is_finite() );
            assert( massmatrix_pseudo.is_finite() );
            assert( massmatrix_volume.is_finite() );
            
            LOG << "assemble Hodge matrices..." << nl;
            
            SparseMatrix hodgematrix_scalar = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
            
            SparseMatrix hodgematrix_vector = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
            
            SparseMatrix hodgematrix_pseudo = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );
            
            SparseMatrix hodgematrix_volume = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r );
            
            assert( hodgematrix_scalar.is_finite() );
            assert( hodgematrix_vector.is_finite() );
            assert( hodgematrix_pseudo.is_finite() );
            assert( hodgematrix_volume.is_finite() );
            
            LOG << "experiments..." << nl;
            
            for( int i = 0; i < number_of_random_trial_samples; i++ )
            {

                const auto& trial_scalarfield = massmatrix_scalar.createinputvector();

                for( int j = 0; j < number_of_random_test_samples; j++ )
                {
                    const auto& test_scalarfield = massmatrix_scalar.createinputvector();

                    // Float mass_errors_scalar[l-l_min][r-r_min] = ;
                    
                    // Float mass_pairing 
                }

                Float mass = 0.;
                
                Assert( mass >= -desired_closeness, mass );
                
                mass_errors_scalar[l-l_min][r-r_min] = std::sqrt( std::abs( mass ) );
                
            }
            
        }
        
        if( l != l_max )
        {
            LOG << "Refinement..." << nl;
        
            M.uniformrefinement();
            
            M.shake_interior_vertices();
        }
        
    } 

    LOG << "Convergence tables" << nl;

    ConvergenceTable contable_scalar;
    ConvergenceTable contable_vector;
    ConvergenceTable contable_pseudo;
    ConvergenceTable contable_volume;
    
    for( int r = r_min; r <= r_max; r++ ) 
    {
        contable_scalar.table_name = "Numerical errors scalar";
        contable_vector.table_name = "Numerical errors vector";
        contable_pseudo.table_name = "Numerical errors pseudo";
        contable_volume.table_name = "Numerical errors volume";

        contable_scalar << printf_into_string("R%d", r-r_min );
        contable_vector << printf_into_string("R%d", r-r_min );
        contable_pseudo << printf_into_string("R%d", r-r_min );
        contable_volume << printf_into_string("R%d", r-r_min );

    }
    
    contable_scalar << nl; 
    contable_vector << nl; 
    contable_pseudo << nl; 
    contable_volume << nl; 

    
    for( int l = l_min; l <= l_max; l++ ) 
    {
        for( int r = r_min; r <= r_max; r++ ) 
        {
            contable_scalar << mass_errors_scalar[l-l_min][r-r_min] / machine_epsilon;
            contable_vector << mass_errors_vector[l-l_min][r-r_min] / machine_epsilon;
            contable_pseudo << mass_errors_pseudo[l-l_min][r-r_min] / machine_epsilon;
            contable_volume << mass_errors_volume[l-l_min][r-r_min] / machine_epsilon;
        }
    
        contable_scalar << nl; 
        contable_vector << nl; 
        contable_pseudo << nl; 
        contable_volume << nl; 
    
    }
        
    contable_scalar.lg(); 
    LOG << "                   " << nl;
    contable_vector.lg(); 
    LOG << "                   " << nl;
    contable_pseudo.lg(); 
    LOG << "                   " << nl;
    contable_volume.lg(); 
    
    
    
    
    
    const Float threshold = 0.1;

    LOG << "Check that differences are below: " << threshold << nl;
    
    for( int l      = l_min; l      <=      l_max; l++      ) 
    for( int r      = r_min; r      <=      r_max; r++      ) 
    {
        if( r < r_max or l < 2 ) 
            continue;
        
        continue; // TODO(martinlicht): find a meaningful test here 
        Assert( mass_errors_scalar[l-l_min][r-r_min] < threshold, mass_errors_scalar[l-l_min][r-r_min], threshold );
        Assert( mass_errors_vector[l-l_min][r-r_min] < threshold, mass_errors_vector[l-l_min][r-r_min], threshold );
        Assert( mass_errors_pseudo[l-l_min][r-r_min] < threshold, mass_errors_pseudo[l-l_min][r-r_min], threshold );
        Assert( mass_errors_volume[l-l_min][r-r_min] < threshold, mass_errors_volume[l-l_min][r-r_min], threshold );
    }
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
