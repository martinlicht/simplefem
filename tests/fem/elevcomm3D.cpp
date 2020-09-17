

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
// #include "../../operators/composed.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
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
        
        auto M = UnitCube3D();
        
        M.check();
        


        const int r_min = 0;
        
        const int r_max = 3;
        
        const int l_min = 0;
        
        const int l_max = 2;
        
        const int number_of_samples = 3;
        
           
        
        Float errors[ M.getinnerdimension()+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        
        
            
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            for( int k = 0; k <= M.getinnerdimension(); k++ ) 
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                cout << "...assemble matrices: l=" << l << " r=" << r << endl;
        
                SparseMatrix elevation_r_1 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r  , 1 );
                SparseMatrix elevation_r_2 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r+1, 1 );
                SparseMatrix elevation_r_3 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r+2, 1 );
                SparseMatrix elevation_r_g = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r  , 3 );
                
                errors[k][ l ][ r ] = 0.;
                
                for( int i = 0; i < number_of_samples; i++ ){

                    auto field = elevation_r_g.createinputvector();
                    field.random();
                    field.normalize();
                    
                    assert( field.isfinite() );
                    
                    const auto path_direct   = elevation_r_g * field;
                    
                    const auto path_indirect = elevation_r_3 * elevation_r_2 * elevation_r_1 * field;
                    
                    const auto error_mass = ( path_direct - path_indirect ).norm();
                    
                    errors[k][l-l_min][r-r_min] = maximum( errors[k][l-l_min][r-r_min], error_mass );
                    
                }
                
            }

            cout << "Refinement..." << endl;
        
            M.uniformrefinement();
            
        } 
        
        
        
        cout << "Convergence tables" << nl;
    
        ConvergenceTable contable[ M.getinnerdimension()+1 ];
        
        for( int k = 0; k <= M.getinnerdimension(); k++ ) 
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
                contable[k] << errors[k][l-l_min][r-r_min];
            
            contable[k] << nl; 
            
        }
        
        
        
        for( int k = 0; k <= M.getinnerdimension(); k++ ) 
        {
            contable[k].print( cout ); 
            cout << "-------------------" << nl;
        }
        
        
        
        cout << "Check that differences are small" << nl;
        
        for( int l      = l_min; l <=                 l_max; l++ ) 
        for( int r      = r_min; r <=                 r_max; r++ ) 
        for( int k      =     0; k <= M.getinnerdimension(); k++ ) 
        {
            assert( errors[k][l-l_min][r-r_min] < 10e-14 );
        }
        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
