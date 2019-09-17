

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../operators/productoperator.hpp"
// #include "../../operators/composed.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.mesh2D.hpp"
#include "../../solver/crm.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.lagrangeincl.hpp"
#include "../../fem/foo.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test for Solution of Neumann Problem" << endl;
        
        cout << std::setprecision(10);

        if(true){

            cout << "Case 2D" << endl;
            
            cout << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = UnitSquare2D();
            
            M.check();
            
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_rhs;
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_sol;


            
            // std::function<FloatVector(const FloatVector&) scalarfield = 
            
            Float xfeq = 1.;
            Float yfeq = 1.;
            

            experiments_sol.push_back( 
                [xfeq,yfeq](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector({ std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ) });
                }
            );

            experiments_rhs.push_back( 
                [xfeq,yfeq](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        xfeq*xfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
                        +
                        yfeq*yfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
                     });
                }
            );

            

            assert( experiments_sol.size() == experiments_rhs.size() );

            cout << "Solving Poisson Problem with Neumann boundary conditions" << endl;
            
            for( int l = 0; l <= 10; l++ ){
                
                for( int r = 1; r <= 1; r++ ) 
                {
                    
                    cout << "...assemble scalar mass matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    SparseMatrix scalar_massmatrix_fac = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 0, r );
                    
                    cout << "...assemble vector mass matrix" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    
                    cout << "...assemble differential matrix and transpose" << endl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    cout << "...assemble inclusion matrix and transpose" << endl;
            
                    SparseMatrix incmatrix = FEECLagrangeInclusionMatrix( M, M.getinnerdimension(), r );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    cout << "...assemble stiffness matrix" << endl;
            
                    // ProductOperator 
                    // auto stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
                    auto op1 = incmatrix_t * diffmatrix_t;
                    auto op2 = op1 * vector_massmatrix;
                    auto op3 = op2 * diffmatrix;
                    auto stiffness = op3 * incmatrix;

                    for( int i = 0; i < experiments_sol.size(); i++){

                        const auto& function_sol = experiments_sol[i];
                        const auto& function_rhs = experiments_rhs[i];
                        
                        cout << "...interpolate explicit solution and rhs" << endl;
            
                        FloatVector interpol_sol = Interpolation( M, M.getinnerdimension(), 0, r, function_sol );
                        FloatVector interpol_rhs = Interpolation( M, M.getinnerdimension(), 0, r, function_rhs );
                        FloatVector interpol_one = Interpolation( M, M.getinnerdimension(), 0, r, constant_one );
                        
                        cout << "...measure kernel component" << endl;
            
                        Float average_sol = interpol_one * ( scalar_massmatrix * interpol_sol );
                        Float average_rhs = interpol_one * ( scalar_massmatrix * interpol_rhs );
                        
                        cout << average_sol << space << average_rhs << endl;

                        cout << "...compute norms of solution and right-hand side:" << endl;
            
                        Float sol_norm = ( scalar_massmatrix_fac * interpol_sol ).norm();
                        Float rhs_norm = ( scalar_massmatrix_fac * interpol_rhs ).norm();
                        
                        cout << "error:    " << sol_norm << endl;
                        cout << "residual: " << rhs_norm << endl;

                        cout << "...create RHS vector" << endl;
            
                        cout << "...calculation (todo)" << endl;
            
                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( M.count_simplices(0), 0. );

                        ConjugateResidualMethod CRM( stiffness );
                        CRM.tolerance = 1e-15;
                        CRM.solve( sol, rhs );

                        cout << "...compute error and residual:" << endl;
            
                        Float errornorm = ( scalar_massmatrix_fac * ( interpol_sol - incmatrix * sol ) ).norm();
                        Float residualnorm = ( rhs - stiffness * sol ).norm();
                        // Float errornorm1 = interpol_sol * ( scalar_massmatrix * interpol_sol );
                        // Float errornorm2 = power( ( scalar_massmatrix_fac * interpol_sol ).norm(), 2. );

                        cout << "error:    " << errornorm << endl;
                        cout << "residual: " << residualnorm << endl;


                        {
                    
                            fstream fs( "./poissonneumann.vtk", std::fstream::out );
                
                            VTK_MeshWriter_Mesh2D vtk( M, fs );
                            vtk.writePreamble( "Poisson-Neumann problem" );
                            vtk.writeCoordinateBlock();
                            vtk.writeTopDimensionalCells();
                            
                            vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                            
                            fs.close();
                    
                        }


                    }
                    
                }

                cout << "Refinement..." << endl;
            
                M.uniformrefinement();
                
                

            } 
        
        }
        
        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
