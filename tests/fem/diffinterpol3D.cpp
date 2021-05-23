#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: (3D) exterior derivative and interpolation";// << endl;
        
        LOG << std::setprecision(10);

        
        LOG << "Initial mesh...";// << endl;
        
        MeshSimplicial3D M = StandardCube3D();
        
        M.check();
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_function;
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_exterior;
        
        experiments_scalar_function.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ std::exp( Constants::pi * vec[0] * vec[1] * vec[2] ) });
            }
        );

        experiments_scalar_exterior.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector( { 
                        Constants::pi * vec[1] * vec[2] * std::exp( Constants::pi * vec[0] * vec[1] * vec[2] ),
                        Constants::pi * vec[0] * vec[2] * std::exp( Constants::pi * vec[0] * vec[1] * vec[2] ),
                        Constants::pi * vec[0] * vec[1] * std::exp( Constants::pi * vec[0] * vec[1] * vec[2] ) 
                    });
            }
        );

        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_vector_function;
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_vector_exterior;
        
        experiments_vector_function.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector( { 
                        -vec[1],
                        vec[0]*vec[1],
                        vec[2] 
                    });
            }
        );

        experiments_vector_exterior.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector( { 
                        0, 
                        vec[0], 
                        1.
                    });
            }
        );


        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_pseudo_function;
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_pseudo_exterior;
        
        experiments_pseudo_function.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector( { 
                        vec[0],
                        vec[1]*vec[1],
                        vec[2]*vec[2]*vec[2] 
                    });
            }
        );

        experiments_pseudo_exterior.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector( { 
                        1. - 2. * vec[1] + 3. * vec[2]*vec[2]
                    });
            }
        );


        
        const int r_min = 1;
        
        const int r_max = 3;
        
        const int l_min = 0;
        
        const int l_max = 2;
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

            
        
        Float errors_scalar[ experiments_scalar_function.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_vector[ experiments_vector_function.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_pseudo[ experiments_pseudo_function.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];

            
        for( int l = l_min; l <= l_max; l++ ){
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                LOG << "...assemble matrices: l=" << l << " r=" << r;// << endl;
        
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                
                SparseMatrix pseudo_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );
                
                SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r-1 );
                
                SparseMatrix scalar_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                SparseMatrix vector_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );

                SparseMatrix pseudo_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 2, r );

                
                for( int i = 0; i < experiments_scalar_function.size(); i++ ){

                    const auto& original_function = experiments_scalar_function[i];
                    const auto& original_exterior = experiments_scalar_exterior[i];
                    
                    FloatVector interpol_function = Interpolation( M, M.getinnerdimension(), 0, r,   original_function );
                    FloatVector interpol_exterior = Interpolation( M, M.getinnerdimension(), 1, r-1, original_exterior );
                    
                    auto commutator_error = interpol_exterior - scalar_diffmatrix * interpol_function;
                    
                    assert( commutator_error.isfinite() );
                    
                    Float commutator_error_mass = commutator_error * ( vector_massmatrix * commutator_error );
                    
                    assert( std::isfinite( commutator_error_mass ) );
                    
                    errors_scalar[i][l-l_min][r-r_min] = std::sqrt( std::fabs( commutator_error_mass ) );
            
                }
                
                for( int i = 0; i < experiments_vector_function.size(); i++ ){

                    const auto& original_function = experiments_vector_function[i];
                    const auto& original_exterior = experiments_vector_exterior[i];
                    
                    FloatVector interpol_function = Interpolation( M, M.getinnerdimension(), 1, r,   original_function );
                    FloatVector interpol_exterior = Interpolation( M, M.getinnerdimension(), 2, r-1, original_exterior );
                    
                    auto commutator_error = interpol_exterior - vector_diffmatrix * interpol_function;
                    
                    assert( commutator_error.isfinite() );
                    
                    Float commutator_error_mass = commutator_error * ( pseudo_massmatrix * commutator_error );
                    
                    assert( std::isfinite( commutator_error_mass ) );
                    
                    errors_vector[i][l-l_min][r-r_min] = std::sqrt( std::fabs( commutator_error_mass ) );
            
                }
                
                for( int i = 0; i < experiments_pseudo_function.size(); i++ ){

                    const auto& original_function = experiments_pseudo_function[i];
                    const auto& original_exterior = experiments_pseudo_exterior[i];
                    
                    FloatVector interpol_function = Interpolation( M, M.getinnerdimension(), 2, r,   original_function );
                    FloatVector interpol_exterior = Interpolation( M, M.getinnerdimension(), 3, r-1, original_exterior );
                    
                    auto commutator_error = interpol_exterior - pseudo_diffmatrix * interpol_function;
                    
                    assert( commutator_error.isfinite() );
                    
                    Float commutator_error_mass = commutator_error * ( volume_massmatrix * commutator_error );
                    
                    assert( std::isfinite( commutator_error_mass ) );
                    
                    errors_pseudo[i][l-l_min][r-r_min] = std::sqrt( std::fabs( commutator_error_mass ) );
            
                }
                
            }

            LOG << "Refinement...";// << endl;
        
            M.uniformrefinement();
            
            

        } 
        
        
        LOG << "Convergence tables" << nl;
    
        ConvergenceTable contable_scalar[ experiments_scalar_function.size() ];
        ConvergenceTable contable_vector[ experiments_vector_function.size() ];
        ConvergenceTable contable_pseudo[ experiments_pseudo_function.size() ];
        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                for( int i = 0; i < experiments_scalar_function.size(); i++ ) 
                    contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min]; // assert( errors_scalar[i][l-l_min][r-r_min] >= 0. ); //
            
                for( int i = 0; i < experiments_vector_function.size(); i++ ) 
                    contable_vector[i] << errors_vector[i][l-l_min][r-r_min]; // assert( errors_vector[i][l-l_min][r-r_min] >= 0. ); //
            
                for( int i = 0; i < experiments_pseudo_function.size(); i++ ) 
                    contable_pseudo[i] << errors_pseudo[i][l-l_min][r-r_min]; // assert( errors_vector[i][l-l_min][r-r_min] >= 0. ); //
            
            }
            
            for( int i = 0; i < experiments_scalar_function.size(); i++ ) contable_scalar[i] << nl; 
            for( int i = 0; i < experiments_vector_function.size(); i++ ) contable_vector[i] << nl; 
            for( int i = 0; i < experiments_pseudo_function.size(); i++ ) contable_pseudo[i] << nl; 
            
        }
            
        for( int i = 0; i < experiments_scalar_function.size(); i++ ) contable_scalar[i].lg(); 
        LOG << "-------------------" << nl;
        for( int i = 0; i < experiments_vector_function.size(); i++ ) contable_vector[i].lg(); 
        LOG << "-------------------" << nl;
        for( int i = 0; i < experiments_pseudo_function.size(); i++ ) contable_pseudo[i].lg(); 
        
        
        
        
        
//         TODO : check for convergence        
//         LOG << "Check that differences are small" << nl;
//         
//         for( int l = l_min; l <= l_max; l++ ) 
//         for( int r = r_min; r <= r_max; r++ ) 
//         {
//             for( int i = 0; i < experiments_scalar_function.size(); i++ ) 
//                 assert( errors_scalar[i][l-l_min][r-r_min] < 10e-6 );
//             
//             for( int i = 0; i < experiments_vector_function.size(); i++ ) 
//                 assert( errors_vector[i][l-l_min][r-r_min] < 10e-6 );
// 
//             for( int i = 0; i < experiments_pseudo_function.size(); i++ ) 
//                 assert( errors_pseudo[i][l-l_min][r-r_min] < 10e-6 );
//         }
        
        
        LOG << "Finished Unit Test";// << endl;
        
        return 0;
}
