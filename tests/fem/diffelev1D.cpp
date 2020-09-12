

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/crm.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.lagrangeincl.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test for Commutativity of Exterior Derivative" << endl;
        
        cout << std::setprecision(10);

        
        cout << "Case 1D" << endl;
        
        cout << "Initial mesh..." << endl;
        
        MeshSimplicial1D M = StandardInterval1D();
        
        M.check();
        
//         assert( M.getouterdimension() == 1 );
        
        cout << "Prepare scalar fields for testing..." << endl;
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> fields;
        
        fields.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ std::exp( vec[0] ) });
            }
        );



        
        

        const int r_min = 1;
        
        const int r_max = 3;
        
        const int l_min = 0;
        
        const int l_max = 4;
        
        
        const int r_plus = 2;
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

            
        
        Float errors[ M.getinnerdimension() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        
            
        for( int l = l_min; l <= l_max; l++ ){
            
            for( int k =     0; k <  M.getinnerdimension(); k++ ) 
            for( int r = r_min; r <= r_max;                 r++ ) 
            {
                
                cout << "...assemble matrices: l=" << l << " k=" << k << " r=" << r << endl;
        
                SparseMatrix lower_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), k, r          );

                SparseMatrix upper_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), k, r + r_plus );

                SparseMatrix diyi_elevation = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k  , r  , r_plus );
                
                SparseMatrix dier_elevation = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k+1, r-1, r_plus );
                
                SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), k+1, r + r_plus - 1 );
                
                
                
                FloatVector interpol_function = Interpolation( M, M.getinnerdimension(), k, r, fields[k] );

                auto path1 = dier_elevation * lower_diffmatrix * interpol_function;

                auto path2 = upper_diffmatrix * diyi_elevation * interpol_function;

                auto commutator_error = path1 - path2;
                
                Float commutator_error_mass = commutator_error * ( massmatrix * commutator_error );

                assert( std::isfinite( commutator_error_mass ) );
                
                errors[k][l-l_min][r-r_min] = std::sqrt( std::fabs( commutator_error_mass ) );
            
                
                
            }

            cout << "Refinement..." << endl;
        
            M.uniformrefinement();
            
            

        } 
        
        
        cout << "Convergence tables" << nl;
    
        ConvergenceTable contables[ M.getinnerdimension() ];
        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                for( int k = 0; k < M.getinnerdimension(); k++ ) {
//                     assert( 7.1 != errors[k][l-l_min][r-r_min] ); 
//                     std::cout << k << space << l << space << r;
                    contables[k] << errors[k][l-l_min][r-r_min];
                }
                        
            }
            
            for( int k = 0; k < M.getinnerdimension(); k++ ) 
                contables[k] << nl; 
            
        }
            
        for( int k = 0; k < M.getinnerdimension(); k++ ) {
            contables[k].print( cout );
            cout << "-------------------" << nl;
        }

        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
