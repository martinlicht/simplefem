#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "(1D) exterior derivative and interpolation" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
        
        LOG << std::setprecision(10);

        
        LOG << "Initial mesh..." << endl;
        
        MeshSimplicial1D M = StandardInterval1D();
        
        M.check();
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_function;
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_exterior;
        
        experiments_scalar_function.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ std::exp( Constants::pi * vec[0] ) });
            }
        );

        experiments_scalar_exterior.push_back( 
            [](const FloatVector& vec) -> FloatVector{
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
            
            LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << endl;
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << endl;

                LOG << "assemble matrices..." << endl;
        
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                
                SparseMatrix scalar_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                
                LOG << "...experiments" << endl;
        
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
                                
            }

            if( l != l_max )
            {
                LOG << "Refinement..." << endl;
            
                M.uniformrefinement();
            }
            
            

        } 
        
        
        LOG << "Convergence tables" << nl;
    
        ConvergenceTable contable_scalar[ experiments_scalar_function.size() ];
        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                for( int i = 0; i < experiments_scalar_function.size(); i++ ) 
                    contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min]; // assert( errors_scalar[i][l-l_min][r-r_min] >= 0. ); //
            
            }
            
            for( int i = 0; i < experiments_scalar_function.size(); i++ ) contable_scalar[i] << nl; 
            
        }
            
        for( int i = 0; i < experiments_scalar_function.size(); i++ ) contable_scalar[i].lg(); 
        
        
        
        
//         TODO : check for convergence        
//         LOG << "Check that differences are small" << nl;
//         
//         for( int l = l_min; l <= l_max; l++ ) 
//         for( int r = r_min; r <= r_max; r++ ) 
//         {
//             for( int i = 0; i < experiments_scalar_function.size(); i++ ) 
//                 assert( errors_scalar[i][l-l_min][r-r_min] < 10e-6 );            
//         }
        
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
