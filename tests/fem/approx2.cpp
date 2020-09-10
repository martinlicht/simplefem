

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


using namespace std;

int main()
{
        
        cout << "Unit Test for Approximating L2 norm of scalar field" << endl;
        
        cout << std::setprecision(10);

        {

            cout << "Case 3D" << endl;
            
            cout << "Initial mesh..." << endl;
            
            MeshSimplicial3D M = UnitCube3D();
            
            M.check();
            
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_field;
            std::vector<Float>                                          experiments_value;


            
            // std::function<FloatVector(const FloatVector&) scalarfield = 
            
            experiments_field.push_back( 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    return FloatVector({ 1.,0.,0. });
                }
            );

            experiments_value.push_back(1.);
            
            experiments_field.push_back( 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    return FloatVector({ 1.,-3.,-2. });
                }
            );

            experiments_value.push_back(14.);
            

            
            experiments_field.push_back( 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    Float x = vec[0]; Float y = vec[1]; Float z = vec[2]; 
                    return FloatVector({ std::sin( x*y ), std::cos( z ), std::exp(x+y+z) });
                }
            );

            experiments_value.push_back(33.42616007376754121867328484399462245509699700639766495311);
            
            
            
            

            
            cout << "Calculate L2 products" << endl;
            
            for( int l = 0; l <= 2; l++ ){
                
                cout << "Refinement..." << endl;
            
                M.uniformrefinement();
                
                cout << "...assemble matrices" << endl;
            
                for( int r = 0; r <= 3; r++ ) 
                {
                    SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );
                    
                    SparseMatrix massmatrix_rhs = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 2, r );
                    
                    for( int i = 0; i < experiments_field.size(); i++){

                        const auto& scalarfield = experiments_field[i];
                        const auto& should_be   = experiments_value[i];

                        FloatVector interpol = Interpolation( M, M.getinnerdimension(), 2, r, scalarfield );

                        Float mass1 = interpol.scalarproductwith( massmatrix * interpol );
                        Float mass2 = power_numerical( ( massmatrix_rhs * interpol ).norm(), 2. );

                        cout << "[ i, l, r ] = [" << i << ", "  << l << ", " << r << "]\t";
                        
                        cout << "mass: " << mass1 << ' ' << mass2 << ' ' << should_be << endl;

                    }
                    
                }
            } 
        
        }
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
