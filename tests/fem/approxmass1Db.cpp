

/**/

#include <cmath>

#include <functional>
#include <vector>
#include <string>

#include "../../base/include.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: (1D) masses are correctly approximated: mass of reference interpolation" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial1D M = UnitInterval1D();
    
    M.check();
    
    
    
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
    
    experiments_scalar_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 1 );
            return FloatVector({ std::exp( vec[0] ) });
        }
    );

    experiments_scalar_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 1 );
            auto x = vec[0];
            return FloatVector({ 
                ( x*x - Float(1.) ) / ( x*x + Float(1.) ) 
            });
        }
    );
    
    
    
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field;
    
    experiments_volume_field = experiments_scalar_field;
    
    
    
    const int r_min = 0;
    
    const int r_max = 5;
    
    const int l_min = 0;
    
    const int l_max = 4;
    
    const int r_ref = 9;
    
    Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    Float errors_volume[ experiments_volume_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
    
    
    
    for( int l = 0; l < l_min; l++ )
        M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ ){
        
        LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
    
        LOG << "...assemble mass matrices" << nl;
    
        SparseMatrix massmatrix_scalar = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r_ref );
        
        SparseMatrix massmatrix_volume = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r_ref );
            
        assert( massmatrix_scalar.is_finite() );
        assert( massmatrix_volume.is_finite() );
            
        for( int r = r_min; r <= r_max; r++ ) 
        {
            LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;

            LOG << "...assemble degree elevation matrices" << nl;
    
            SparseMatrix elevation_scalar = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_ref - r );
            
            SparseMatrix elevation_volume = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r, r_ref - r );
            
            assert( elevation_scalar.is_finite() );
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
            
            for( int i = 0; i < experiments_volume_field.size(); i++ ){

                const auto& volumefield = experiments_volume_field[i];
                
                auto interpol     = Interpolation( M, M.getinnerdimension(), 1, r, volumefield );

                auto interpol_ref = Interpolation( M, M.getinnerdimension(), 1, r_ref, volumefield );
                
                auto error = interpol_ref - elevation_volume * interpol;

                auto error_mass = error * ( massmatrix_volume * error );
                
                Assert( error_mass >= -desired_closeness, error_mass );
                
                errors_volume[i][l-l_min][r-r_min] = std::sqrt( std::abs( error_mass ) );
                
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
    ConvergenceTable contable_volume[ experiments_volume_field.size() ];
    
    for( int r = r_min; r <= r_max; r++ ) 
    {
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
            contable_scalar[i].table_name = "Numerical errors scalar E" + std::to_string(i);
        for( int i = 0; i < experiments_volume_field.size(); i++ ) 
            contable_volume[i].table_name = "Numerical errors volume E" + std::to_string(i);

        for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
            contable_scalar[i] << printf_into_string("R%d", r-r_min );
        for( int i = 0; i < experiments_volume_field.size(); i++ ) 
            contable_volume[i] << printf_into_string("R%d", r-r_min );

    }
    for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
    for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 

    
    for( int l = l_min; l <= l_max; l++ ) 
    {
        
        for( int r = r_min; r <= r_max; r++ ) 
        {
            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min];
        
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                contable_volume[i] << errors_volume[i][l-l_min][r-r_min];
        
        }
        
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
        for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
        
    }
        
    for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i].lg(); 
    LOG << "                   " << nl;
    for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i].lg(); 
    
    
    
    
    
    LOG << "Check that differences are below: " << desired_closeness_for_sqrt << nl;
    
    for( int l      = l_min; l      <=      l_max; l++      ) 
    for( int r      = r_min; r      <=      r_max; r++      ) 
    {
        if( r < r_max or l < 3 ) 
            continue;
        
        // continue; // TODO(martinlicht): find a meaningful test here 
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
            Assert( errors_scalar[i][l-l_min][r-r_min] < desired_closeness_for_sqrt, errors_scalar[i][l-l_min][r-r_min], desired_closeness_for_sqrt );
        
        for( int i = 0; i < experiments_volume_field.size(); i++ )
            Assert( errors_volume[i][l-l_min][r-r_min] < desired_closeness_for_sqrt, errors_volume[i][l-l_min][r-r_min], desired_closeness_for_sqrt );
    }
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
