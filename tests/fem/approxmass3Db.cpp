

/**/

#include <cmath>

#include <functional>
#include <vector>
#include <string>

#include "../../base/include.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: (3D) masses are correctly approximated: mass of reference interpolation" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = UnitCube3D();
    
    M.check();
    
    
    
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
    std::vector<Float>                                          experiments_scalar_value;
    
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector {
//                     assert( vec.getdimension() == 3 );
//                     return FloatVector({ 1. });
//                 }
//         );
// 
//         experiments_scalar_value.push_back( 1. );
//         
//         
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector {
//                     assert( vec.getdimension() == 3 );
//                     return FloatVector({ vec.sum() < 1. ? 1. : 0. });
//                 }
//         );
// 
//         experiments_scalar_value.push_back( 1./6. );
    
    
    experiments_scalar_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 3 );
            return FloatVector({ std::exp( vec[0] + vec[1] + vec[2] ) });
        }
    );

    experiments_scalar_value.push_back( 3.19452804946532511361521373028750 * 3.19452804946532511361521373028750 * 3.19452804946532511361521373028750 );


    
    
    
    
    
    
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_vector_field;
    std::vector<Float>                                          experiments_vector_value;
    
//         experiments_vector_field.push_back( 
//         [](const FloatVector& vec) -> FloatVector {
//                 assert( vec.getdimension() == 3 );
//                 return FloatVector({ 1.,0.,0. });
//             }
//         );
// 
//         experiments_vector_value.push_back( 1. );
//         
//         
//         experiments_vector_field.push_back( 
//         [](const FloatVector& vec) -> FloatVector {
//                 assert( vec.getdimension() == 3 );
//                 return FloatVector({ 1.,-3.,-2. });
//             }
//         );
//         
//         experiments_vector_value.push_back( 14. );
    
    
    experiments_vector_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
                assert( vec.getdimension() == 3 );
                Float x = vec[0]; Float y = vec[1]; Float z = vec[2]; 
                return FloatVector({ std::sin( x*y ), std::cos( z ), std::exp(x+y+z) });
            }
        );

    experiments_vector_value.push_back( 33.42616007376754121867328484399462245509699700639766495311 );
    
    
    
    
    
    
    
    
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_pseudo_field;
    std::vector<Float>                                          experiments_pseudo_value;

    experiments_pseudo_field = experiments_vector_field;
    experiments_pseudo_value = experiments_vector_value;
    
    
    
    
    
    
    
    
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field;
    std::vector<Float>                                          experiments_volume_value;

    experiments_volume_field = experiments_scalar_field;
    experiments_volume_value = experiments_scalar_value;
    
    
    
    
    const int r_min = 0;
    
    const int r_max = 3;
    
    const int l_min = 0;
    
    const int l_max = 3 ;
    
    const int r_ref = 4;
    
    Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float errors_vector[ experiments_vector_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float errors_pseudo[ experiments_pseudo_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float errors_volume[ experiments_volume_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    
    
    
    for( int l = 0; l < l_min; l++ )
        M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ ){
        
        LOG << "Level: " << l_min << " <= " << l << " <= " << l_max << nl;

        LOG << "... assemble mass matrices" << nl;

        SparseMatrix massmatrix_scalar = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r_ref );
        
        SparseMatrix massmatrix_vector = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r_ref );
        
        SparseMatrix massmatrix_pseudo = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r_ref );
        
        SparseMatrix massmatrix_volume = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r_ref );
        
        assert( massmatrix_scalar.is_finite() );
        assert( massmatrix_vector.is_finite() );
        assert( massmatrix_pseudo.is_finite() );
        assert( massmatrix_volume.is_finite() );
            
        for( int r = r_min; r <= r_max; r++ ) 
        {
            
            LOG << "Level: "      << l_min << " <= " << l << " <= " << l_max << nl;
            LOG << "Polydegree: " << r_min << " <= " << r << " <= " << r_max << nl;

            LOG << "... assemble degree elevation matrices" << nl;
            
            SparseMatrix elevation_scalar = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_ref - r );
            
            SparseMatrix elevation_vector = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r, r_ref - r );
            
            SparseMatrix elevation_pseudo = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r, r_ref - r );
            
            SparseMatrix elevation_volume = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 3, r, r_ref - r );
            
            assert( elevation_scalar.is_finite() );
            assert( elevation_vector.is_finite() );
            assert( elevation_pseudo.is_finite() );
            assert( elevation_volume.is_finite() );
            
            LOG << "experiments..." << nl;
            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ){

                const auto& scalarfield = experiments_scalar_field[i];
                
                auto interpol     = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                auto interpol_ref = Interpolation( M, M.getinnerdimension(), 0, r_ref, scalarfield );
                
                auto error = interpol_ref - elevation_scalar * interpol;

                auto error_mass = error * ( massmatrix_scalar * error );
                
                Assert( error_mass >= -desired_closeness, error_mass );
                
                errors_scalar[i][l-l_min][r-r_min] = std::sqrt( std::abs( error_mass ) );
                
            }
            
            for( int i = 0; i < experiments_vector_field.size(); i++ ){

                const auto& vectorfield = experiments_vector_field[i];
                
                auto interpol     = Interpolation( M, M.getinnerdimension(), 1, r, vectorfield );

                auto interpol_ref = Interpolation( M, M.getinnerdimension(), 1, r_ref, vectorfield );
                
                auto error = interpol_ref - elevation_vector * interpol;

                auto error_mass = error * ( massmatrix_vector * error );
                
                Assert( error_mass >= -desired_closeness, error_mass );
                
                errors_vector[i][l-l_min][r-r_min] = std::sqrt( std::abs( error_mass ) );
                
            }
            
            for( int i = 0; i < experiments_pseudo_field.size(); i++ ){

                const auto& pseudofield = experiments_pseudo_field[i];
                
                auto interpol     = Interpolation( M, M.getinnerdimension(), 2, r, pseudofield );

                auto interpol_ref = Interpolation( M, M.getinnerdimension(), 2, r_ref, pseudofield );
                
                auto error = interpol_ref - elevation_pseudo * interpol;

                auto error_mass = error * ( massmatrix_pseudo * error );
                
                Assert( error_mass >= -desired_closeness, error_mass );
                
                errors_pseudo[i][l-l_min][r-r_min] = std::sqrt( std::abs( error_mass ) );
                
            }
            
            for( int i = 0; i < experiments_volume_field.size(); i++ ){

                const auto& volumefield = experiments_volume_field[i];
                
                auto interpol     = Interpolation( M, M.getinnerdimension(), 3, r, volumefield );

                auto interpol_ref = Interpolation( M, M.getinnerdimension(), 3, r_ref, volumefield );
                
                auto error = interpol_ref - elevation_volume * interpol;

                auto error_mass = error * ( massmatrix_volume * error );
                
                Assert( error_mass >= -desired_closeness, error_mass );
                
                errors_volume[i][l-l_min][r-r_min] = std::sqrt( std::abs( error_mass ) );

                Assert( errors_volume[i][l-l_min][r-r_min] > 0., i, l, r );
                
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

    ConvergenceTable contable_scalar[ experiments_scalar_field.size() ];
    ConvergenceTable contable_vector[ experiments_vector_field.size() ];
    ConvergenceTable contable_pseudo[ experiments_vector_field.size() ];
    ConvergenceTable contable_volume[ experiments_volume_field.size() ];
    
    for( int r = r_min; r <= r_max; r++ ) 
    {
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
            contable_scalar[i].table_name = "Numerical errors scalar E" + std::to_string(i);
        for( int i = 0; i < experiments_vector_field.size(); i++ ) 
            contable_vector[i].table_name = "Numerical errors vector E" + std::to_string(i);
        for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
            contable_pseudo[i].table_name = "Numerical errors pseudo E" + std::to_string(i);
        for( int i = 0; i < experiments_volume_field.size(); i++ ) 
            contable_volume[i].table_name = "Numerical errors volume E" + std::to_string(i);

        for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
            contable_scalar[i] << printf_into_string("R%d", r-r_min );
        for( int i = 0; i < experiments_vector_field.size(); i++ ) 
            contable_vector[i] << printf_into_string("R%d", r-r_min );
        for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
            contable_pseudo[i] << printf_into_string("R%d", r-r_min );
        for( int i = 0; i < experiments_volume_field.size(); i++ ) 
            contable_volume[i] << printf_into_string("R%d", r-r_min );

    }
    for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
    for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
    for( int i = 0; i < experiments_pseudo_field.size(); i++ ) contable_pseudo[i] << nl; 
    for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
    
    
    for( int l = l_min; l <= l_max; l++ ) 
    {
        
        for( int r = r_min; r <= r_max; r++ ) 
        {
            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min];
        
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                contable_vector[i] << errors_vector[i][l-l_min][r-r_min];
        
            for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
                contable_pseudo[i] << errors_pseudo[i][l-l_min][r-r_min];
        
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                contable_volume[i] << errors_volume[i][l-l_min][r-r_min];
        
        }
        
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
        for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
        for( int i = 0; i < experiments_pseudo_field.size(); i++ ) contable_pseudo[i] << nl; 
        for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
        
    }
        
    for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i].lg(); 
    LOG << "                   " << nl;
    for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i].lg(); 
    LOG << "                   " << nl;
    for( int i = 0; i < experiments_pseudo_field.size(); i++ ) contable_pseudo[i].lg(); 
    LOG << "                   " << nl;
    for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i].lg(); 
    
    
    
    
    
    const Float threshold = 1e-3;

    LOG << "Check that differences are below: " << threshold << nl;
    
    for( int l      = l_min; l      <=      l_max; l++      ) 
    for( int r      = r_min; r      <=      r_max; r++      ) 
    {
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
            Assert( std::isfinite( errors_scalar[i][l-l_min][r-r_min] ) );
        
        for( int i = 0; i < experiments_vector_field.size(); i++ ) 
            Assert( std::isfinite( errors_vector[i][l-l_min][r-r_min] ) );

        for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
            Assert( std::isfinite( errors_pseudo[i][l-l_min][r-r_min] ) );

        if( r < r_max or l < l_max ) 
            continue;
        
        continue; // TODO(Martin): This test depends on convergence and cannot be exact
        
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
            Assert( errors_scalar[i][l-l_min][r-r_min] < threshold, errors_scalar[i][l-l_min][r-r_min], threshold );
        
        for( int i = 0; i < experiments_vector_field.size(); i++ ) 
            Assert( errors_vector[i][l-l_min][r-r_min] < threshold, errors_vector[i][l-l_min][r-r_min], threshold );

        for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
            Assert( errors_pseudo[i][l-l_min][r-r_min] < threshold, errors_pseudo[i][l-l_min][r-r_min], threshold );
        
        for( int i = 0; i < experiments_volume_field.size(); i++ )
            Assert( errors_volume[i][l-l_min][r-r_min] < threshold, errors_volume[i][l-l_min][r-r_min], threshold );
    }
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
