

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
        
        cout << "Unit Test for Approximating L2 norms of fields" << endl;
        
        cout << std::setprecision(10);

        
        
        
        
        
        cout << "Initial mesh..." << endl;
        
        MeshSimplicial1D M = UnitInterval1D();
        
        M.check();
        
        cout << "Prepare scalar fields for testing..." << endl;
        

        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
        
        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 1 );
                return FloatVector({ std::exp( vec[0] ) });
            }
        );

        

        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field;
        
        experiments_volume_field = experiments_scalar_field;
        



        
        
        const int r_min = 0;
        
        const int r_max = 5;
        
        const int l_min = 0;
        
        const int l_max = 10;
        
        
//         const int r_plus = 2;
        
        const int R = 5;
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        
        
        Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_volume[ experiments_volume_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        
        
        
        
        cout << "Calculate L2 products" << endl;
        
        for( int l = l_min; l <= l_max; l++ ){
            
            cout << "Refinement... l = " << l << endl;
        
            M.uniformrefinement();
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                cout << "...assemble matrices" << endl;
        
                SparseMatrix massmatrix_scalar_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, R );
                
                SparseMatrix massmatrix_volume_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, R );
                
                SparseMatrix elevation_scalar = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, R - r );

                SparseMatrix elevation_volume = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r, R - r  );

                
                for( int i = 0; i < experiments_scalar_field.size(); i++){

                    const auto& scalarfield = experiments_scalar_field[i];
        
                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                    FloatVector interpol_plus = Interpolation( M, M.getinnerdimension(), 0, R, scalarfield );
                    
                    FloatVector error = interpol_plus - elevation_scalar * interpol;

                    Float error_mass = error * ( massmatrix_scalar_plus * error );
                    
                    errors_scalar[i][l-l_min][r-r_min] = std::sqrt( error_mass );
                    
                }
                
                for( int i = 0; i < experiments_volume_field.size(); i++){

                    const auto& volumefield = experiments_volume_field[i];
        
                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 1, r, volumefield );

                    FloatVector interpol_plus = Interpolation( M, M.getinnerdimension(), 1, R, volumefield );
                    
                    FloatVector error = interpol_plus - elevation_volume * interpol;

                    Float error_mass = error * ( massmatrix_volume_plus * error );
                    
                    errors_volume[i][l-l_min][r-r_min] = std::sqrt( error_mass );
                    
                }
                
            }
        } 
    
        cout << "Convergence tables" << nl;
    
        ConvergenceTable contable_scalar[ experiments_scalar_field.size() ];
        ConvergenceTable contable_volume[ experiments_volume_field.size() ];
        
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
            
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i].print( cout ); 
        cout << "-------------------" << nl;
        for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i].print( cout ); 
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
