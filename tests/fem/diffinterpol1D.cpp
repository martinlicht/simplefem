
#include <cmath>

#include <functional>
#include <vector>
#include <string>

#include "../../base/include.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: (1D) exterior derivative and interpolation" << nl;
    
    
    MeshSimplicial1D M = StandardInterval1D();
    
    M.check();
    
    
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_function;
    std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_exterior;
    
    experiments_scalar_function.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            return FloatVector({ std::exp( Constants::pi * vec[0] ) });
        }
    );

    experiments_scalar_exterior.push_back( 
        [](const FloatVector& vec) -> FloatVector {
            assert( vec.getdimension() == 2 );
            return FloatVector( { 
                    Constants::pi * std::exp( Constants::pi * vec[0] ),
                    0., 
                });
        }
    );
    
    
    
    const int r_min = 1;
    
    const int r_max = 4;
    
    const int l_min = 0;
    
    const int l_max = 8;
    
    
    for( int l = 0; l < l_min; l++ )
        M.uniformrefinement();

        
    
    Float errors_scalar[ experiments_scalar_function.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        
    for( int l = l_min; l <= l_max; l++ ){
        
        LOG << "Level: " << l_min << " <= " << l << " <= " << l_max << nl;
        
        for( int r = r_min; r <= r_max; r++ ) 
        {
            
            LOG << "Level: "      << l_min << " <= " << l << " <= " << l_max << nl;
            LOG << "Polydegree: " << r_min << " <= " << r << " <= " << r_max << nl;

            LOG << "assemble matrices..." << nl;
    
            SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
            
            SparseMatrix scalar_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

            
            LOG << "... experiments" << nl;
    
            for( int i = 0; i < experiments_scalar_function.size(); i++ ){

                const auto& original_function = experiments_scalar_function[i];
                const auto& original_exterior = experiments_scalar_exterior[i];
                
                FloatVector interpol_function = Interpolation( M, M.getinnerdimension(), 0, r,   original_function );
                FloatVector interpol_exterior = Interpolation( M, M.getinnerdimension(), 1, r-1, original_exterior );
                
                auto commutator_error = interpol_exterior - scalar_diffmatrix * interpol_function;
                
                assert( commutator_error.is_finite() );
                
                Float commutator_error_mass = commutator_error * ( vector_massmatrix * commutator_error );
                
                assert( std::isfinite( commutator_error_mass ) );
                Assert( commutator_error_mass >= -desired_closeness, commutator_error_mass );
                
                errors_scalar[i][l-l_min][r-r_min] = std::sqrt( std::fabs( commutator_error_mass ) );
        
            }
                            
        }

        if( l != l_max )
        {
            LOG << "Refinement..." << nl;
        
            M.uniformrefinement();

            // M.shake_interior_vertices(); // The inner and outer dimension may differ.
        }
        
        

    } 
    
    
    LOG << "Convergence tables" << nl;

    ConvergenceTable contable_scalar[ experiments_scalar_function.size() ];
    
    for( int r = r_min; r <= r_max; r++ ) 
    {
        for( int i = 0; i < experiments_scalar_function.size(); i++ ) 
            contable_scalar[i].table_name = "Numerical errors scalar E" + std::to_string(i);
        
        for( int i = 0; i < experiments_scalar_function.size(); i++ ) 
            contable_scalar[i] << printf_into_string("R%d", r );

    }
    for( int i = 0; i < experiments_scalar_function.size(); i++ ) contable_scalar[i] << nl; 

    
    
    for( int l = l_min; l <= l_max; l++ ) 
    {
        
        for( int r = r_min; r <= r_max; r++ ) 
        {
            
            for( int i = 0; i < experiments_scalar_function.size(); i++ ) 
                contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min]; // Assert( errors_scalar[i][l-l_min][r-r_min] >= -desired_closeness ); //
        
        }
        
        for( int i = 0; i < experiments_scalar_function.size(); i++ ) contable_scalar[i] << nl; 
        
    }
        
    for( int i = 0; i < experiments_scalar_function.size(); i++ ) contable_scalar[i].lg(); 
    
    
    
    
    LOG << "Check that differences are below: " << desired_closeness_for_sqrt << nl;
    
    for( int l = l_min; l <= l_max; l++ ) 
    for( int r = r_min; r <= r_max; r++ ) 
    {
        if( r < r_max or l < l_max ) continue;

        // TODO(Martin): This test depends on convergence and cannot be exact

        for( int i = 0; i < experiments_scalar_function.size(); i++ ) 
            Assert( errors_scalar[i][l-l_min][r-r_min] < desired_closeness_for_sqrt, errors_scalar[i][l-l_min][r-r_min], desired_closeness_for_sqrt, r, l );
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
