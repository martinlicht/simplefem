

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
// #include "../../operators/composed.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/crm.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.lagrangeincl.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test: (1D) degree elevations commute" << endl;
        
        cout << std::setprecision(10);

        
        cout << "Initial mesh..." << endl;
        
        MeshSimplicial1D M = UnitInterval1D();
        
        M.check();
        
        


        const int r_min = 0;
        
        const int r_max = 5;
        
        const int l_min = 0;
        
        const int l_max = 2;
        
        
        const int number_of_trials = 3;
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

            
        
        Float errors_scalar[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_volume[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        
        
            
        for( int l = l_min; l <= l_max; l++ ){
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                cout << "...assemble matrices: l=" << l << " r=" << r << endl;
        
                SparseMatrix scalar_elevation_r_1 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r  , 1 );
                SparseMatrix scalar_elevation_r_2 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r+1, 1 );
                SparseMatrix scalar_elevation_r_3 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r+2, 1 );
                SparseMatrix scalar_elevation_r_g = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r  , 3 );
                
                SparseMatrix volume_elevation_r_1 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r  , 1 );
                SparseMatrix volume_elevation_r_2 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r+1, 1 );
                SparseMatrix volume_elevation_r_3 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r+2, 1 );
                SparseMatrix volume_elevation_r_g = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r  , 3 );
                
                errors_scalar[ l ][ r ] = 0.;
                errors_volume[ l ][ r ] = 0.;
        
                for( int i = 0; i < number_of_trials; i++){

                    auto scalarfield = scalar_elevation_r_g.createinputvector();
                    scalarfield.random();
                    
                    const auto path_direct   = scalar_elevation_r_g * scalarfield;
                    
                    const auto path_indirect = scalar_elevation_r_3 * scalar_elevation_r_2 * scalar_elevation_r_1 * scalarfield;
                    
                    const auto error_mass = ( path_direct - path_indirect ).norm();
                    
                    errors_scalar[l-l_min][r-r_min] = maximum( errors_scalar[l-l_min][r-r_min], error_mass );
                    
                }
                
                for( int i = 0; i < number_of_trials; i++){

                    auto volumefield = volume_elevation_r_g.createinputvector();
                    volumefield.random();
                    
                    const auto path_direct   = volume_elevation_r_g * volumefield;
                    
                    const auto path_indirect = volume_elevation_r_3 * volume_elevation_r_2 * volume_elevation_r_1 * volumefield;
                    
                    const auto error_mass = ( path_direct - path_indirect ).norm();
                    
                    errors_volume[l-l_min][r-r_min] = maximum( errors_volume[l-l_min][r-r_min], error_mass );
                    
                }
                
                
                
            }

            cout << "Refinement..." << endl;
        
            M.uniformrefinement();
            
        } 
        
        
        cout << "Convergence tables" << nl;
    
        ConvergenceTable contable_scalar;
        ConvergenceTable contable_volume;
        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                contable_scalar << errors_scalar[l-l_min][r-r_min];
            
                contable_volume << errors_volume[l-l_min][r-r_min];
            
            }
            
            contable_scalar << nl; 
            contable_volume << nl; 
            
        }
        
        
        
        cout << "Convergence tables: scalars" << nl;
        {
            contable_scalar.print( cout ); 
            cout << "-------------------" << nl;
        }
        
        cout << "Convergence tables: volumes" << nl;
        {
            contable_volume.print( cout ); 
            cout << "-------------------" << nl;
        }
        
        
        
        
        
        cout << "Check that differences are small" << nl;
        
        for( int l      = l_min; l      <=      l_max; l++      ) 
        for( int r      = r_min; r      <=      r_max; r++      ) 
        {
            assert( errors_scalar[l-l_min][r-r_min] < 10e-14 );
            assert( errors_volume[l-l_min][r-r_min] < 10e-14 );
        }

        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
