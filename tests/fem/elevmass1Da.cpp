

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test: (1D) degree elevation of interpolation has the mass of higher order interpolation" << endl;
        
        cout << std::setprecision(10);

        
        
        
        
        
        cout << "Initial mesh..." << endl;
        
        MeshSimplicial1D M = StandardInterval1D();
        
        M.check();
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
        
        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ std::exp( vec[0] ) });
            }
        );

        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field;
        
        experiments_volume_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ std::exp( vec[0] ), 0. });
            }
        );

        
        
        


        
        
        const int r_min = 0;
        
        const int r_max = 5;
        
        const int l_min = 0;
        
        const int l_max = 4;
        
        
        const int r_plus_max = 3;
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        
        
        Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        Float errors_volume[ experiments_volume_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        
        
        
        
        cout << "Calculate L2 products" << endl;
        
        for( int l = l_min; l <= l_max; l++ ){
            
            cout << "Numerical calculations..." << endl;
            
            for( int r      = r_min; r      <=      r_max; r++      ) 
            for( int r_plus =     0; r_plus <= r_plus_max; r_plus++ ) 
            {
                
                SparseMatrix massmatrix_scalar = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                
                SparseMatrix massmatrix_volume = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
                
                SparseMatrix massmatrix_scalar_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r + r_plus );
                
                SparseMatrix massmatrix_volume_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r + r_plus );
                
                SparseMatrix elevation_scalar = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_plus );

                SparseMatrix elevation_volume = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r, r_plus );

                
                for( int i = 0; i < experiments_scalar_field.size(); i++){

                    const auto& scalarfield = experiments_scalar_field[i];
        
                    FloatVector interpol      = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                    FloatVector interpol_elev = elevation_scalar * interpol;
                    
                    Float mass      = interpol * ( massmatrix_scalar * interpol );

                    Float mass_elev = interpol_elev * ( massmatrix_scalar_plus * interpol_elev );

                    Float error_mass = mass - mass_elev;
                    
                    errors_scalar[i][l-l_min][r-r_min][r_plus] = error_mass;
                    
                }
                
                for( int i = 0; i < experiments_volume_field.size(); i++){

                    const auto& volumefield = experiments_volume_field[i];
        
                    FloatVector interpol      = Interpolation( M, M.getinnerdimension(), 1, r, volumefield );

                    FloatVector interpol_elev = elevation_volume * interpol;
                    
                    Float mass      = interpol * ( massmatrix_volume * interpol );

                    Float mass_elev = interpol_elev * ( massmatrix_volume_plus * interpol_elev );

                    Float error_mass = mass - mass_elev;
                    
                    errors_volume[i][l-l_min][r-r_min][r_plus] = error_mass;
                    
                }
                
            }
            
            cout << "Refinement..." << endl;
        
            M.uniformrefinement();
            
        } 
    
        cout << "Convergence tables" << nl;
    
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
        
        
        
        cout << "Convergence tables: scalars" << nl;
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
        {
            contable_scalar[i].print( cout ); 
            cout << "-------------------" << nl;
        }
        
        cout << "Convergence tables: volumes" << nl;
        for( int i = 0; i < experiments_volume_field.size(); i++ )
        {
            contable_volume[i].print( cout ); 
            cout << "-------------------" << nl;
        }
        
        
        
        
        
        cout << "Check that differences are small" << nl;
        
        for( int l      = l_min; l      <=      l_max; l++      ) 
        for( int r      = r_min; r      <=      r_max; r++      ) 
        for( int r_plus =     0; r_plus <= r_plus_max; r_plus++ ) 
        {
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                assert( errors_scalar[i][l-l_min][r-r_min][r_plus] < 10e-14 );
            
            for( int i = 0; i < experiments_volume_field.size(); i++ )
                assert( errors_volume[i][l-l_min][r-r_min][r_plus] < 10e-14 );
        }
            
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
