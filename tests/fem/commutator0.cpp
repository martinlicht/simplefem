

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
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.lagrangeincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test for Approximation of Gradients" << endl;
        
        cout << std::setprecision(10);

        {
            MeshSimplicial2D M = UnitTriangle2D();

            SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, 1 );

            std::cout << signpower(0) << space << signpower(1) << space << signpower(2) << space << signpower(-1) << space << diffmatrix << endl;

            // return 0;
        }


        if(true){

            cout << "Case 2D" << endl;
            
            cout << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_func;
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_grad;
            

            
            Float xfeq = 1.;

            experiments_func.push_back( 
                [&](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ ( vec[0] > 0 and vec[1] > 0 ) ? xfeq * vec[0] : 0. });
                }
            );

            experiments_grad.push_back( 
                [&](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector( { 
                            ( vec[0] > 0 and vec[1] > 0 ) ? 1. : 0. ,
                            ( vec[0] > 0 and vec[1] > 0 ) ? 0. : 0. ,
                        });
                }
            );

            // experiments_sol.push_back( 
            //     [xfeq,yfeq](const FloatVector& vec) -> FloatVector{
            //         assert( vec.getdimension() == 2 );
            //         return FloatVector({ std::cos( xfeq * Constants::pi * vec[0] ) });
            //     }
            // );

            // experiments_grad.push_back( 
            //     [xfeq,yfeq](const FloatVector& vec) -> FloatVector{
            //         assert( vec.getdimension() == 2 );
            //         return FloatVector( { 
            //                 -xfeq * Constants::pi * std::sin( xfeq * Constants::pi * vec[0] ),
            //                 0., 
            //             });
            //     }
            // );


            cout << "Solving Poisson Problem with Neumann boundary conditions" << endl;
            
            for( int l = 0; l <= 10; l++ ){
                
                for( int r = 1; r <= 3; r++ ) 
                {
                    
                    cout << "...assemble scalar mass matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    SparseMatrix scalar_massmatrix_fac = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 0, r );
                    
                    cout << "...assemble vector mass matrix" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    
                    SparseMatrix vector_massmatrix_fac = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 1, r-1 );
                    
                    cout << "...assemble differential matrix and transpose" << endl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    // cout << "...assemble inclusion matrix and transpose" << endl;
            
                    // SparseMatrix incmatrix = LagrangeInclusionMatrix( M, M.getinnerdimension(), r );

                    // SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    for( int i = 0; i < experiments_func.size(); i++){

                        const auto& function_func = experiments_func[i];
                        const auto& function_grad = experiments_grad[i];
                        
                        cout << "...interpolate function and gradient" << endl;
            
                        FloatVector interpol_func = Interpolation( M, M.getinnerdimension(), 0, r,   function_func );
                        FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r-1, function_grad );
                        
                        cout << "...measure interpolation commutativity" << endl;
            
                        Float commutatorerror = ( vector_massmatrix_fac * ( interpol_grad - ( diffmatrix * interpol_func ) ) ).norm();
                        cout << "commutator error: " << commutatorerror << endl;
                        
                        
                        
                        
                        // {
                    
                        //     fstream fs( "./poissonneumann.vtk", std::fstream::out );
                
                        //     VTKWriter vtk( M, fs );
                        //     vtk.writePreamble( "Poisson-Neumann problem" );
                        //     vtk.writeCoordinateBlock();
                        //     vtk.writeTopDimensionalCells();
                            
                        //     vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                        //     vtk.writeCellVectorData( interpol_grad, "gradient_interpolation" , 0.1 );
                            
                        //     fs.close();
                    
                        // }


                    }
                    
                }

                cout << "Refinement..." << endl;
            
                M.uniformrefinement();
                
                

            } 
        
        }
        
        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
