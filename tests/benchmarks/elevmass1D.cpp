

/**/

#include "../base/include.hpp"
#include "../dense/densematrix.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/mesh.simplicial1D.hpp"
#include "../mesh/examples1D.hpp"
#include "../fem/global.massmatrix.hpp"
#include "../fem/global.elevation.hpp"
#include "../fem/utilities.hpp"
#include "../utility/convergencetable.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: (1D) degree elevation of interpolation has the mass of higher order interpolation" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial1D M = StandardInterval1D();
    
    M.check();
    
    
    
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
    
    experiments_scalar_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            return FloatVector({ std::exp( vec[0] ) });
        }
    );

    experiments_scalar_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            auto x = vec[0];
            return FloatVector({ 
                ( x*x - 1. ) / ( x*x + 1. ) 
            });
        }
    );

    
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field;
    
    experiments_volume_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            return FloatVector({ std::exp( vec[0] ), 0. });
        }
    );
    
    experiments_volume_field.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            auto x = vec[0];
            return FloatVector({ 
                ( x*x - 1. ) / ( x*x + 1. ) ,
                0.
            });
        }
    );

    
    
    
    const int r_min = 0;
    
    const int r_max = 4;
    
    const int l_min = 0;
    
    const int l_max = 4;
    
    const int r_plus_max = 3;
    
    Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
    Float errors_volume[ experiments_volume_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
    
    
    
    for( int l = 0; l < l_min; l++ )
        M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ ){
        
        LOG << "Numerical calculations..." << nl;
        
        for( int r      = r_min; r      <=      r_max; r++      ) 
        for( int r_plus =     0; r_plus <= r_plus_max; r_plus++ ) 
        {
            
            LOG << "Level: "             << l_min << " <= " << l << " <= " << l_max << nl;
            LOG << "Polynomial degree: " << r_min << " <= " << r << " <= " << r_max << nl;
            LOG << "additional degree: " <<   "0" << " <= " << r_plus << " <= " << r_plus_max << nl;
            
            SparseMatrix massmatrix_scalar_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r + r_plus );
            
            SparseMatrix massmatrix_volume_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r + r_plus );
            
            SparseMatrix elevation_scalar = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_plus );

            SparseMatrix elevation_volume = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r, r_plus );

            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ){

                const auto& scalarfield = experiments_scalar_field[i];
    
                FloatVector interpol = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                FloatVector interpol_plus = Interpolation( M, M.getinnerdimension(), 0, r + r_plus, scalarfield );
                
                FloatVector error = interpol_plus - elevation_scalar * interpol;

                Float error_mass = error * ( massmatrix_scalar_plus * error );
                
                errors_scalar[i][l-l_min][r-r_min][r_plus] = std::sqrt( error_mass );
                
            }
            
            for( int i = 0; i < experiments_volume_field.size(); i++ ){

                const auto& volumefield = experiments_volume_field[i];
    
                FloatVector interpol = Interpolation( M, M.getinnerdimension(), 1, r, volumefield );

                FloatVector interpol_plus = Interpolation( M, M.getinnerdimension(), 1, r + r_plus, volumefield );
                
                FloatVector error = interpol_plus - elevation_volume * interpol;

                Float error_mass = error * ( massmatrix_volume_plus * error );
                
                errors_volume[i][l-l_min][r-r_min][r_plus] = std::sqrt( error_mass );
                
            }
            
        }
        
        LOG << "Refinement..." << nl;
    
        M.uniformrefinement();
        
    } 

    LOG << "Convergence tables for the case of largest degree jump" << nl;

    ConvergenceTable contable_scalar[ experiments_scalar_field.size() ];
    ConvergenceTable contable_volume[ experiments_volume_field.size() ];
    
    for( int l = l_min; l <= l_max; l++ ) 
    {
        
        for( int r = r_min; r <= r_max; r++ ) 
        {
            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min][r_plus_max];
        
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                contable_volume[i] << errors_volume[i][l-l_min][r-r_min][r_plus_max];
        
        }
        
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
        for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
        
    }
    
    
    
    LOG << "Convergence tables: scalars" << nl;
    for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
    {
        LOG << contable_scalar[i].text(); 
        LOG << "-------------------" << nl;
    }
    
    LOG << "Convergence tables: volumes" << nl;
    for( int i = 0; i < experiments_volume_field.size(); i++ )
    {
        LOG << contable_volume[i].text(); 
        LOG << "-------------------" << nl;
    }
    
    
    
    
    
    
//         LOG << "Check that differences are small: " << desired_closeness << nl;
//         
//         for( int l      = l_min; l      <=      l_max; l++      ) 
//         for( int r      = r_min; r      <=      r_max; r++      ) 
//         for( int r_plus =     0; r_plus <= r_plus_max; r_plus++ ) 
//         {
//             for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
//                 assert( errors_scalar[i][l-l_min][r-r_min][r_plus] < desired_closeness );
//             
//             
//             for( int i = 0; i < experiments_volume_field.size(); i++ )
//                 assert( errors_volume[i][l-l_min][r-r_min][r_plus] < desired_closeness );
//         }
        
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
