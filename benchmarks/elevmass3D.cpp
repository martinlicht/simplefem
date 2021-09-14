

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../basic.hpp"
#include "../dense/densematrix.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/mesh.simplicial3D.hpp"
#include "../mesh/examples3D.hpp"
#include "../fem/local.polynomialmassmatrix.hpp"
#include "../fem/global.massmatrix.hpp"
#include "../fem/global.elevation.hpp"
#include "../fem/utilities.hpp"
#include "../utility/convergencetable.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test: (3D) degree elevation of interpolation has the mass of higher order interpolation" << endl;
        
        LOG << std::setprecision(10);

        LOG << "Initial mesh..." << endl;
        
        MeshSimplicial3D M = StandardCube3D();
        
        M.check();
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
        std::vector<Float>                                          experiments_scalar_value;
        
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                     assert( vec.getdimension() == 3 );
//                     return FloatVector({ 1. });
//                 }
//         );
// 
//         experiments_scalar_value.push_back( 1. );
//         
//         
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                     assert( vec.getdimension() == 3 );
//                     return FloatVector({ vec.sum() < 1. ? 1. : 0. });
//                 }
//         );
// 
//         experiments_scalar_value.push_back( 1./6. );
        
        
        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ std::exp( vec[0] + vec[1] + vec[2] ) });
            }
        );

        experiments_scalar_value.push_back( 3.19452804946532511361521373028750 * 3.19452804946532511361521373028750 * 3.19452804946532511361521373028750 );


        
        
        
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_vector_field;
        std::vector<Float>                                          experiments_vector_value;
        
//         experiments_vector_field.push_back( 
//         [](const FloatVector& vec) -> FloatVector{
//                 assert( vec.getdimension() == 3 );
//                 return FloatVector({ 1.,0.,0. });
//             }
//         );
// 
//         experiments_vector_value.push_back( 1. );
//         
//         
//         experiments_vector_field.push_back( 
//         [](const FloatVector& vec) -> FloatVector{
//                 assert( vec.getdimension() == 3 );
//                 return FloatVector({ 1.,-3.,-2. });
//             }
//         );
//         
//         experiments_vector_value.push_back( 14. );
        
        
        experiments_vector_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
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
        
        const int r_max = 2;
        
        const int l_min = 0;
        
        const int l_max = 2;
        
        const int r_plus_max = 3;
        
        Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        Float errors_vector[ experiments_vector_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        Float errors_pseudo[ experiments_pseudo_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        Float errors_volume[ experiments_volume_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            LOG << "Numerical calculations..." << endl;
            
            for( int r      = r_min; r      <=      r_max; r++      ) 
            for( int r_plus =     0; r_plus <= r_plus_max; r_plus++ ) 
            {
                
                SparseMatrix massmatrix_scalar_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r + r_plus );
                
                SparseMatrix massmatrix_vector_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r + r_plus );
                
                SparseMatrix massmatrix_pseudo_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r + r_plus );
                
                SparseMatrix massmatrix_volume_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r + r_plus );
                
                SparseMatrix elevation_scalar       = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_plus );

                SparseMatrix elevation_vector       = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r, r_plus );

                SparseMatrix elevation_pseudo       = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r, r_plus );

                SparseMatrix elevation_volume       = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 3, r, r_plus );

                
                for( int i = 0; i < experiments_scalar_field.size(); i++ ){

                    const auto& scalarfield = experiments_scalar_field[i];
        
                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                    FloatVector interpol_plus = Interpolation( M, M.getinnerdimension(), 0, r + r_plus, scalarfield );
                    
                    FloatVector error = interpol_plus - elevation_scalar * interpol;

                    Float error_mass = error * ( massmatrix_scalar_plus * error );
                    
                    errors_scalar[i][l-l_min][r-r_min][r_plus] = std::sqrt( error_mass );
                    
                }
                
                for( int i = 0; i < experiments_vector_field.size(); i++ ){

                    const auto& vectorfield = experiments_vector_field[i];
        
                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 1, r, vectorfield );

                    FloatVector interpol_plus = Interpolation( M, M.getinnerdimension(), 1, r + r_plus, vectorfield );
                    
                    FloatVector error = interpol_plus - elevation_vector * interpol;

                    Float error_mass = error * ( massmatrix_vector_plus * error );
                    
                    errors_vector[i][l-l_min][r-r_min][r_plus] = std::sqrt( error_mass );
                    
                }
                
                for( int i = 0; i < experiments_pseudo_field.size(); i++ ){

                    const auto& pseudofield = experiments_pseudo_field[i];
        
                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 2, r, pseudofield );

                    FloatVector interpol_plus = Interpolation( M, M.getinnerdimension(), 2, r + r_plus, pseudofield );
                    
                    FloatVector error = interpol_plus - elevation_pseudo * interpol;

                    Float error_mass = error * ( massmatrix_pseudo_plus * error );
                    
                    errors_pseudo[i][l-l_min][r-r_min][r_plus] = std::sqrt( error_mass );
                    
                }
                
                for( int i = 0; i < experiments_volume_field.size(); i++ ){

                    const auto& volumefield = experiments_volume_field[i];
        
                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 3, r, volumefield );

                    FloatVector interpol_plus = Interpolation( M, M.getinnerdimension(), 3, r + r_plus, volumefield );
                    
                    FloatVector error = interpol_plus - elevation_volume * interpol;

                    Float error_mass = error * ( massmatrix_volume_plus * error );
                    
                    errors_volume[i][l-l_min][r-r_min][r_plus] = std::sqrt( error_mass );
                    
                }
                
            }
            
            LOG << "Refinement..." << endl;
        
            M.uniformrefinement();
            
        } 
    
        LOG << "Convergence tables for the case of largest degree jump" << nl;
    
        ConvergenceTable contable_scalar[ experiments_scalar_field.size() ];
        ConvergenceTable contable_vector[ experiments_vector_field.size() ];
        ConvergenceTable contable_pseudo[ experiments_vector_field.size() ];
        ConvergenceTable contable_volume[ experiments_volume_field.size() ];
        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                    contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min][r_plus_max];
            
                for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                    contable_vector[i] << errors_vector[i][l-l_min][r-r_min][r_plus_max];
            
                for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
                    contable_pseudo[i] << errors_pseudo[i][l-l_min][r-r_min][r_plus_max];
            
                for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                    contable_volume[i] << errors_volume[i][l-l_min][r-r_min][r_plus_max];
            
            }
            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
            for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
            for( int i = 0; i < experiments_pseudo_field.size(); i++ ) contable_pseudo[i] << nl; 
            for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
            
        }
        
        
        
        LOG << "Convergence tables: scalars" << nl;
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
        {
            LOG << contable_scalar[i].text(); 
            LOG << "-------------------" << nl;
        }
        
        LOG << "Convergence tables: vectors" << nl;
        for( int i = 0; i < experiments_vector_field.size(); i++ ) 
        {
            LOG << contable_vector[i].text(); 
            LOG << "-------------------" << nl;
        }
        
        LOG << "Convergence tables: pseudos" << nl;
        for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
        {
            LOG << contable_pseudo[i].text(); 
            LOG << "-------------------" << nl;
        }
        
        LOG << "Convergence tables: volumes" << nl;
        for( int i = 0; i < experiments_volume_field.size(); i++ )
        {
            LOG << contable_volume[i].text(); 
            LOG << "-------------------" << nl;
        }
        
        
        
        
        
        
//         LOG << "Check that differences are small" << nl;
//         
//         for( int l      = l_min; l      <=      l_max; l++      ) 
//         for( int r      = r_min; r      <=      r_max; r++      ) 
//         for( int r_plus =     0; r_plus <= r_plus_max; r_plus++ ) 
//         {
//             for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
//                 assert( errors_scalar[i][l-l_min][r-r_min][r_plus] < 10e-14 );
//             
//             for( int i = 0; i < experiments_vector_field.size(); i++ ) 
//                 assert( errors_vector[i][l-l_min][r-r_min][r_plus] < 10e-14 );
//             
//             for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
//                 assert( errors_pseudo[i][l-l_min][r-r_min][r_plus] < 10e-14 );
//             
//             for( int i = 0; i < experiments_volume_field.size(); i++ )
//                 assert( errors_volume[i][l-l_min][r-r_min][r_plus] < 10e-14 );
//         }
            
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
