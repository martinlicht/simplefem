

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
// #include "../../operators/composed.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples2D.hpp"
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
        
        cout << "Unit Test for Commutativity of Exterior Derivative" << endl;
        
        cout << std::setprecision(10);

        
        cout << "Case 2D" << endl;
        
        cout << "Initial mesh..." << endl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        cout << "Prepare scalar fields for testing..." << endl;
        


        const int r_min = 0;
        
        const int r_max = 5;
        
        const int l_min = 0;
        
        const int l_max = 2;
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

            
        
        Float errors_scalar[ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_vector[ l_max - l_min + 1 ][ r_max - r_min + 1 ];

            
        for( int l = l_min; l <= l_max; l++ ){
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                cout << "...assemble matrices: l=" << l << " r=" << r << endl;
        
                SparseMatrix scalar_elevation_r_1 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r  , 1 );
                SparseMatrix scalar_elevation_r_2 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r+1, 1 );
                SparseMatrix scalar_elevation_r_3 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r+2, 1 );
                SparseMatrix scalar_elevation_r_g = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r  , 3 );
                
                SparseMatrix vector_elevation_r_1 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r  , 1 );
                SparseMatrix vector_elevation_r_2 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r+1, 1 );
                SparseMatrix vector_elevation_r_3 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r+2, 1 );
                SparseMatrix vector_elevation_r_g = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r  , 3 );
                
                SparseMatrix volume_elevation_r_1 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r  , 1 );
                SparseMatrix volume_elevation_r_2 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r+1, 1 );
                SparseMatrix volume_elevation_r_3 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r+2, 1 );
                SparseMatrix volume_elevation_r_g = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r  , 3 );
                
                
                
                
            }

            cout << "Refinement..." << endl;
        
            M.uniformrefinement();
            
            

        } 
        
        
//         cout << "Convergence tables" << nl;
//     
//         ConvergenceTable contable_scalar;
//         ConvergenceTable contable_vector;
//         ConvergenceTable contable_vector;
//         
//         for( int l = l_min; l <= l_max; l++ ) 
//         {
//             
//             for( int r = r_min; r <= r_max; r++ ) 
//             {
//                 
//                 for( int i = 0; i < experiments_scalar_function.size(); i++ ) 
//                     contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min];
//             
//                 for( int i = 0; i < experiments_vector_function.size(); i++ ) 
//                     contable_vector[i] << errors_vector[i][l-l_min][r-r_min];
//             
//             }
//             
//             for( int i = 0; i < experiments_scalar_function.size(); i++ ) contable_scalar[i] << nl; 
//             for( int i = 0; i < experiments_vector_function.size(); i++ ) contable_vector[i] << nl; 
//             
//         }
//             
//         for( int i = 0; i < experiments_scalar_function.size(); i++ ) contable_scalar[i].print( cout ); 
//         cout << "-------------------" << nl;
//         for( int i = 0; i < experiments_vector_function.size(); i++ ) contable_vector[i].print( cout ); 

        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
