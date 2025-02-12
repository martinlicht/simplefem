

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test for Approximating L2 norms of fields" << endl;
        
        cout << std::setprecision(10);

        
        
        
        
        
        cout << "Initial mesh..." << endl;
        
        MeshSimplicial3D M = UnitCube3D();
        
        M.check();
        
        cout << "Prepare scalar fields for testing..." << endl;
        

        
        
        
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
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        
        Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_vector[ experiments_vector_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_pseudo[ experiments_pseudo_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_volume[ experiments_volume_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        
        
        
        cout << "Calculate L2 products" << endl;
        
        for( int l = l_min; l <= l_max; l++ ){
            
            cout << "Numerical calculations..." << endl;
        
            M.uniformrefinement();
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                cout << "...assemble matrices" << endl;
        
                SparseMatrix massmatrix_scalar = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                
                SparseMatrix massmatrix_vector = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
                
                SparseMatrix massmatrix_pseudo = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );
                
                SparseMatrix massmatrix_volume = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r );
                
                assert( massmatrix_scalar.isfinite() );
                assert( massmatrix_vector.isfinite() );
                assert( massmatrix_pseudo.isfinite() );
                assert( massmatrix_volume.isfinite() );
                
                for( int i = 0; i < experiments_scalar_field.size(); i++ ){

                    const auto& scalarfield = experiments_scalar_field[i];
                    const auto& should_be   = experiments_scalar_value[i];

                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                    assert( interpol.isfinite() );

                    Float mass = interpol * ( massmatrix_scalar * interpol );
                    
                    errors_scalar[i][l-l_min][r-r_min] = std::sqrt( std::abs( mass - should_be ) );
                    
                    cout << "mass: " << mass << ' ' << should_be << endl; 

                }
                
                for( int i = 0; i < experiments_vector_field.size(); i++ ){

                    const auto& vectorfield = experiments_vector_field[i];
                    const auto& should_be   = experiments_vector_value[i];

                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 1, r, vectorfield );

                    assert( interpol.isfinite() );

                    Float mass = interpol * ( massmatrix_vector * interpol );
                    
                    errors_vector[i][l][r] = std::sqrt( std::abs( mass - should_be ) );
                    
                    cout << "mass: " << mass << ' ' << should_be << endl; 

                }
                
                for( int i = 0; i < experiments_vector_field.size(); i++ ){

                    const auto& vectorfield = experiments_vector_field[i];
                    const auto& should_be   = experiments_vector_value[i];

                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 2, r, vectorfield );

                    assert( interpol.isfinite() );

                    Float mass = interpol * ( massmatrix_pseudo * interpol );
                    
                    errors_pseudo[i][l][r] = std::sqrt( std::abs( mass - should_be ) );
                    
                    cout << "mass: " << mass << ' ' << should_be << endl; 

                }
                
                for( int i = 0; i < experiments_volume_field.size(); i++ ){

                    const auto& volumefield = experiments_volume_field[i];
                    const auto& should_be   = experiments_volume_value[i];

                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 3, r, volumefield );

                    assert( interpol.isfinite() );

                    Float mass = interpol * ( massmatrix_volume * interpol );
                    
                    errors_volume[i][l][r] = std::sqrt( std::abs( mass - should_be ) );
                    
                    cout << "mass: " << mass << ' ' << should_be << endl; 

                }
                
            }
        } 
    
    
        cout << "Convergence tables" << nl;
    
        ConvergenceTable contable_scalar[ experiments_scalar_field.size() ];
        ConvergenceTable contable_vector[ experiments_vector_field.size() ];
        ConvergenceTable contable_pseudo[ experiments_vector_field.size() ];
        ConvergenceTable contable_volume[ experiments_volume_field.size() ];
        
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
            
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i].print( cout ); 
        cout << "-------------------" << nl;
        for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i].print( cout ); 
        cout << "-------------------" << nl;
        for( int i = 0; i < experiments_pseudo_field.size(); i++ ) contable_pseudo[i].print( cout ); 
        cout << "-------------------" << nl;
        for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i].print( cout ); 
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
