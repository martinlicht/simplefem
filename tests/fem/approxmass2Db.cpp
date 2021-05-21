

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: (2D) masses are correctly approximated: mass of reference interpolation";// << endl;
        
        LOG << std::setprecision(10);

        LOG << "Initial mesh...";// << endl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
        std::vector<Float>                                          experiments_scalar_value;
        
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 assert( vec.getdimension() == 2 );
//                 return FloatVector({ ( vec[0] > 0 and vec[1] > 0 ) ? 1. : 0. });
//             }
//         );
// 
//         experiments_scalar_value.push_back( 1. );
//         
//         
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 assert( vec.getdimension() == 2 );
//                 return FloatVector({ 1. });
//             }
//         );
// 
//         experiments_scalar_value.push_back( 4. );
//         
//         
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 assert( vec.getdimension() == 2 );
//                 return FloatVector({ ( vec[0] * vec[1] > 0 ) ? 1. : 0. });
//             }
//         );
// 
//         experiments_scalar_value.push_back( 2. );
        
        
        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ std::exp( vec[0] + vec[1] ) });
            }
        );

        experiments_scalar_value.push_back( 3.62686040784701876 * 3.62686040784701876 );

        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ ( vec[0] > 0 and vec[1] > 0 ) ? std::exp( vec[0] ) : 0. });
            }
        );

        experiments_scalar_value.push_back(3.194528049465325113615213730287503906590157785275923662043);
        
        
        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ ( vec[0] * vec[1] > 0 ) ? std::exp( vec[0] ) : 0. });
            }
        );

        experiments_scalar_value.push_back(3.626860407847018767668213982801261704886342012321135721309);
        
        
        
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_vector_field;
        std::vector<Float>                                          experiments_vector_value;
        
//         experiments_vector_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 assert( vec.getdimension() == 2 );
//                 // return FloatVector({ ( vec[0]*vec[1] > 0 ) ? 1. : 0., ( vec[0]*vec[1] > 0 ) ? 10. : 0. });
//                 return FloatVector({ 1., 0. });
//             }
//         );
// 
//         experiments_vector_value.push_back( 4. );
        
        
        experiments_vector_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                Float x = vec[0]; Float y = vec[1];
                return FloatVector({ std::exp(x-y), std::exp(x-y) });
            }
        );

        // int_(-1)^(1) int_(-1)^(1) exp(x-y)^2 + exp(x-y)^2 dx dy
        experiments_vector_value.push_back( 26.308232836016486629201989612067059822501324553083772160298096942 );
        

        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field;
        std::vector<Float>                                          experiments_volume_value;

//         experiments_volume_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 assert( vec.getdimension() == 2 );
//                 return FloatVector({ ( vec[0] > 0 and vec[1] > 0 ) ? 1. : 0. });
//             }
//         );
// 
//         experiments_volume_value.push_back( 1. );
//         
//         
//         experiments_volume_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 assert( vec.getdimension() == 2 );
//                 return FloatVector({ ( vec[0] * vec[1] > 0 ) ? 1. : 0. });
//             }
//         );
// 
//         experiments_volume_value.push_back( 2. );
        
        
        experiments_volume_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ ( vec[0] > 0 and vec[1] > 0 ) ? std::exp( vec[0] ) : 0. });
            }
        );

        experiments_volume_value.push_back(3.194528049465325113615213730287503906590157785275923662043);
        
        
        experiments_volume_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ ( vec[0] * vec[1] > 0 ) ? std::exp( vec[0] ) : 0. });
            }
        );

        experiments_volume_value.push_back(3.626860407847018767668213982801261704886342012321135721309);
        
        
        
        
        const int r_min = 0;
        
        const int r_max = 5;
        
        const int l_min = 0;
        
        const int l_max = 4;
        
        const int r_ref = 9;
        
        Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_vector[ experiments_vector_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_volume[ experiments_volume_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            LOG << "Refinement...";// << endl;
        
            M.uniformrefinement();
            
            SparseMatrix massmatrix_scalar = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r_ref );
            
            SparseMatrix massmatrix_vector = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r_ref );
                
            SparseMatrix massmatrix_volume = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r_ref );
                
            assert( massmatrix_scalar.isfinite() );
            assert( massmatrix_vector.isfinite() );
            assert( massmatrix_volume.isfinite() );
                
            for( int r = r_min; r <= r_max; r++ ) 
            {
                LOG << "...assemble matrices";// << endl;
        
                SparseMatrix elevation_scalar = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_ref - r );
                
                SparseMatrix elevation_vector = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r, r_ref - r );
                
                SparseMatrix elevation_volume = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r, r_ref - r );
                
                assert( elevation_scalar.isfinite() );
                assert( elevation_volume.isfinite() );
                
                for( int i = 0; i < experiments_scalar_field.size(); i++ ){

                    const auto& scalarfield = experiments_scalar_field[i];
                    
                    auto interpol     = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                    auto interpol_ref = Interpolation( M, M.getinnerdimension(), 0, r_ref, scalarfield );
                    
                    auto error = interpol_ref - elevation_scalar * interpol;

                    auto error_mass = error * ( massmatrix_scalar * error );
                    
                    errors_scalar[i][l-l_min][r-r_min] = std::sqrt( std::abs( error_mass ) );
                    
                }
                
                for( int i = 0; i < experiments_vector_field.size(); i++ ){

                    const auto& vectorfield = experiments_vector_field[i];
                    
                    auto interpol     = Interpolation( M, M.getinnerdimension(), 1, r, vectorfield );

                    auto interpol_ref = Interpolation( M, M.getinnerdimension(), 1, r_ref, vectorfield );
                    
                    auto error = interpol_ref - elevation_vector * interpol;

                    auto error_mass = error * ( massmatrix_vector * error );
                    
                    errors_vector[i][l-l_min][r-r_min] = std::sqrt( std::abs( error_mass ) );
                    
                }
                
                for( int i = 0; i < experiments_volume_field.size(); i++ ){

                    const auto& volumefield = experiments_volume_field[i];
                    
                    auto interpol     = Interpolation( M, M.getinnerdimension(), 2, r, volumefield );

                    auto interpol_ref = Interpolation( M, M.getinnerdimension(), 2, r_ref, volumefield );
                    
                    auto error = interpol_ref - elevation_volume * interpol;

                    auto error_mass = error * ( massmatrix_volume * error );
                    
                    errors_volume[i][l-l_min][r-r_min] = std::sqrt( std::abs( error_mass ) );
                    
                }
                
            }
        } 
    
        LOG << "Convergence tables" << nl;
    
        ConvergenceTable contable_scalar[ experiments_scalar_field.size() ];
        ConvergenceTable contable_vector[ experiments_vector_field.size() ];
        ConvergenceTable contable_volume[ experiments_volume_field.size() ];
        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                    contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min];
            
                for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                    contable_vector[i] << errors_vector[i][l-l_min][r-r_min];
            
                for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                    contable_volume[i] << errors_volume[i][l-l_min][r-r_min];
            
            }
            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
            for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
            for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
            
        }
            
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i].lg(); 
        LOG << "-------------------" << nl;
        for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i].lg(); 
        LOG << "-------------------" << nl;
        for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i].lg(); 
        
        
        
        
        
        LOG << "Check that differences are small" << nl;
        
        for( int l      = l_min; l      <=      l_max; l++      ) 
        for( int r      = r_min; r      <=      r_max; r++      ) 
        {
            if( r < r_max or l < 3 ) 
                continue;
            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                assert( errors_scalar[i][l-l_min][r-r_min] < 10e-6 );
            
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                assert( errors_vector[i][l-l_min][r-r_min] < 10e-6 );

            for( int i = 0; i < experiments_volume_field.size(); i++ )
                assert( errors_volume[i][l-l_min][r-r_min] < 10e-6 );
        }
        
        
        LOG << "Finished Unit Test";// << endl;
        
        return 0;
}
