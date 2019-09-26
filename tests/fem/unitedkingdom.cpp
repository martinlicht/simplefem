

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../operators/productoperator.hpp"
// #include "../../operators/composed.hpp"
#include "../../dense/densematrix.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.mesh2D.hpp"
#include "../../solver/crm.hpp"
#include "../../solver/pcrm.hpp"
#include "../../solver/resdes.hpp"
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
            
            MeshSimplicial2D M = UnitedKingdom(); //UnitSquare2D();
            
            M.check();
            
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_rhs;
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_grad;
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_sol;

            
            // std::function<FloatVector(const FloatVector&) scalarfield = 
            
            experiments_sol.push_back( 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector({ 
                        vec[0] * vec[0] * vec[1] * vec[1]
                        });
                }
            );

            experiments_grad.push_back( 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector( { 
                            2. * vec[0] * vec[1] * vec[1],
                            2. * vec[1] * vec[0] * vec[0], 
                        });
                }
            );

            experiments_rhs.push_back( 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        2. * vec[0] * vec[0] + 2. * vec[1] * vec[1] 
                     });
                }
            );

            

            cout << "Solving Poisson Problem with Neumann boundary conditions over the United Kingdom" << endl;

            int max_l = 2;
            int max_r = 1;
            
            for( int l = 0; l <= max_l; l++ ){
                
                cout << "Level: " << l << std::endl;
                cout << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                for( int r = 1; r <= max_r; r++ ) 
                {
                    
                    cout << "...assemble scalar mass matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    SparseMatrix scalar_massmatrix_fac = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 0, r );
                    
                    cout << "...assemble vector mass matrix" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    
                    SparseMatrix vector_massmatrix_fac = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 1, r-1 );
                    
                    SparseMatrix vector_massmatrix_fac_t = vector_massmatrix_fac.getTranspose();
                    cout << "...assemble differential matrix and transpose" << endl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    cout << "...assemble inclusion matrix and transpose" << endl;
            
                    SparseMatrix incmatrix = FEECLagrangeInclusionMatrix( M, M.getinnerdimension(), r );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    cout << "...assemble stiffness matrix" << endl;
            
                    // ProductOperator 
                    // auto stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
                    auto op1_t = incmatrix_t * diffmatrix_t;
                    auto op2_t = op1_t * vector_massmatrix_fac_t;
                    auto op1   = diffmatrix * incmatrix;
                    auto op2   = vector_massmatrix_fac * op1;
                    auto stiffness = op2_t * op2;

                    // auto opr1 = diffmatrix & incmatrix;
                    // auto opr  = vector_massmatrix_fac & opr1;
                    // auto opl  = opr.getTranspose(); 
                    // auto stiffness = opl & opr;
                    
                    // stiffness.sortentries();
                    // auto stiffness_csr = MatrixCSR( stiffness );
                    
                    //auto stiffness_invprecon = DiagonalOperator( stiffness.getdimin(), 1. );
                    // auto stiffness_invprecon = InverseDiagonalPreconditioner( stiffness );
                    // std::cout << "Average value of diagonal preconditioner: " << stiffness_invprecon.getdiagonal().average() << std::endl;

                    for( int i = 0; i < experiments_rhs.size(); i++){

                        const auto& function_sol = experiments_sol[i];
                        const auto& function_grad= experiments_grad[i];
                        const auto& function_rhs = experiments_rhs[i];
                        
                        cout << "...interpolate rhs" << endl;
            
                        FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, r,   function_sol  );
                        FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r-1, function_grad );
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        
                        FloatVector interpol_rhs_cells  = Interpolation( M, M.getinnerdimension(), 0, 0,   function_rhs  );
                        
                        FloatVector interpol_one  = Interpolation( M, M.getinnerdimension(), 0, r, constant_one );
                        
                        cout << "...measure kernel component: " << std::flush;
            
                        Float average_sol = interpol_one * ( scalar_massmatrix * interpol_sol );
                        Float average_rhs = interpol_one * ( scalar_massmatrix * interpol_rhs );
                        Float domain_area = interpol_one * ( scalar_massmatrix * interpol_one );

                        cout << average_rhs << endl;

                        interpol_rhs = interpol_rhs - ( average_rhs / domain_area ) * interpol_one;

                        cout << "...measure interpolation commutativity" << endl;
            
                        cout << "...compute norms of solution and right-hand side:" << endl;
            
                        Float rhs_norm = ( scalar_massmatrix_fac * interpol_rhs ).norm();
                        
                        cout << "rhs norm:      " << rhs_norm << endl;

                        cout << "...create RHS vector" << endl;

                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( M.count_simplices(0), 0. );
                        
                        cout << "...iterative solver" << endl;
                        
                        sol.zero();
                        FloatVector eins( sol.getdimension(), 1. );
                        
                        if(false)
                        for( int t = 0; t < 4; t++ )
                        {
                            timestamp start = gettimestamp();
                            ConjugateResidualMethod CRM( stiffness );
                            CRM.print_modulo = 1;//+sol.getdimension()/1000;
                            CRM.tolerance = 1e-15;
                            CRM.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t " << end - start << std::endl;
                            sol = sol - ( interpol_one * ( scalar_massmatrix * ( incmatrix * sol ) ) ) * eins;
                        }

                        for( int t = 0; t < 10; t++ )
                        {
                            timestamp start = gettimestamp();
                            ResidualDescent( stiffness, rhs, sol, sol.getdimension(), 1e-15 );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t " << end - start << std::endl;
                            sol = sol - ( interpol_one * ( scalar_massmatrix * ( incmatrix * sol ) ) ) * eins;
                        }


                        cout << "...compute error and residual:" << endl;
            
                        // Float errornorm     = ( scalar_massmatrix_fac * ( interpol_sol  - incmatrix * sol ) ).norm();
                        // Float graderrornorm = ( vector_massmatrix_fac * ( interpol_grad - diffmatrix * incmatrix * sol ) ).norm();
                        Float residualnorm  = ( rhs - stiffness * sol ).norm();
                        
                        // FloatVector gradfoo = diffmatrix * ( interpol_sol - incmatrix * sol );
                        // Float graderrornorm = gradfoo.scalarproductwith( vector_massmatrix * gradfoo );
                        // Float errornorm1 = interpol_sol * ( scalar_massmatrix * interpol_sol );
                        // Float errornorm2 = power( ( scalar_massmatrix_fac * interpol_sol ).norm(), 2. );

                        // cout << "error:     " << errornorm     << endl;
                        // cout << "graderror: " << graderrornorm << endl;
                        cout << "residual:  " << residualnorm  << endl;

                        {
                    
                            fstream fs( "./poissonneumann.uk.vtk", std::fstream::out );
                
                            VTK_MeshWriter_Mesh2D vtk( M, fs );
                            vtk.writePreamble( "Poisson-Neumann problem" );
                            vtk.writeCoordinateBlock();
                            vtk.writeTopDimensionalCells();
                            
                            vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                            vtk.writeCellScalarData( interpol_rhs_cells, "rhs_interpolation" , 0.1 );
                            
                            fs.close();
                    
                        }


                    }
                    
                }

                cout << "Refinement..." << endl;
            
                if( l != max_l ) M.uniformrefinement();
                
                

            } 
        
        }
        
        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
