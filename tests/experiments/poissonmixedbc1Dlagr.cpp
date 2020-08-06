

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../operators/composedoperators.hpp"
// #include "../../operators/composed.hpp"
#include "../../dense/densematrix.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../vtk/vtkwriter.mesh1D.hpp"
#include "../../solver/minres.hpp"
#include "../../solver/crm.hpp"
#include "../../solver/minres.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.lagrangeincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test for Solution of Neumann Problem" << endl;
        
        cout << std::setprecision(10);

        if(true){

            cout << "Initial mesh..." << endl;
            
            MeshSimplicial1D M = StandardSquare1D();
            
            M.check();
            
            M.set_flag( 0, 0, SimplexFlagDirichlet );
            
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 1 );
                        return FloatVector({ 1. });
                    };
            
            
            
            


            
            // std::function<FloatVector(const std::function<FloatVector(const FloatVector&) ) >scalarfield = 
            
            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 1 );
                    return FloatVector({ 
                        std::exp( vec[0] ) - 1. - Constants::euler * vec[0]
                    });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_grad = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 1 );
                    return FloatVector( { 
                        std::exp( vec[0] )      - Constants::euler 
                    });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 1 );
                    return FloatVector({ 
                        - std::exp( vec[0] ) 
                     });
                };
            

            

            

            cout << "Solving Poisson Problem with Dirichlet boundary conditions" << endl;

            int min_l = 3; 
            int max_l = 15;
            
            ConvergenceTable contable;
            

            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                cout << "Level: " << l << std::endl;
                cout << "# E/V: " << M.count_edges() << "/" << M.count_vertices() << nl;
                
                const int r = 1;
                
                {
                    
                    cout << "...assemble scalar mass matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    cout << "...assemble vector mass matrix" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    
                    cout << "...assemble differential matrix and transpose" << endl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    cout << "...assemble inclusion matrix and transpose" << endl;
            
                    SparseMatrix incmatrix = LagrangeInclusionMatrix( M, M.getinnerdimension(), r );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    cout << "...assemble stiffness matrix" << endl;
            
                    // ProductOperator 
                    // auto stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
                    // auto op1 = incmatrix_t * diffmatrix_t;
                    // auto op2 = op1 * vector_massmatrix;
                    // auto op3 = op2 * diffmatrix;
                    // auto stiffness = op3 * incmatrix;

//                     auto opr1 = diffmatrix & incmatrix;
//                     auto opr  = vector_massmatrix_fac & opr1;
//                     auto opl  = opr.getTranspose(); 
//                     auto stiffness = opl & opr;

                    auto opr  = diffmatrix & incmatrix;
                    auto opl  = opr.getTranspose(); 
                    auto stiffness = opl & ( vector_massmatrix & opr );
                    
                    stiffness.sortentries();
                    auto stiffness_csr = MatrixCSR( stiffness );
                    
                    auto stiffness_invprecon = DiagonalOperator( stiffness.getdimin(), 1. );
                    //auto stiffness_invprecon = InverseDiagonalPreconditioner( stiffness );
                    std::cout << "Average value of diagonal preconditioner: " << stiffness_invprecon.getdiagonal().average() << std::endl;

                    {

                        const auto& function_sol  = experiment_sol;
                        const auto& function_grad = experiment_grad;
                        const auto& function_rhs  = experiment_rhs;
                        
                        cout << "...interpolate explicit solution and rhs" << endl;
            
                        FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, r,   function_sol  );
                        FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r-1, function_grad );
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        
                        cout << "...measure interpolation commutativity" << endl;
            
                        FloatVector commutator = interpol_grad - diffmatrix * interpol_sol;
                        Float commutatorerror = std::sqrt( commutator * ( vector_massmatrix * commutator ) );
                        cout << "commutator error: " << commutatorerror << endl;
                        
                        cout << "...compute norms of solution and right-hand side:" << endl;
            
                        Float sol_norm = std::sqrt( interpol_sol * ( scalar_massmatrix * interpol_sol ) );
                        Float rhs_norm = std::sqrt( interpol_rhs * ( scalar_massmatrix * interpol_rhs ) );
                        
                        cout << "solution norm: " << sol_norm << endl;
                        cout << "rhs norm:      " << rhs_norm << endl;

                        cout << "...create RHS vector" << endl;

                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( M.count_simplices(0), 0. );
                        
                        cout << "...iterative solver" << endl;
                        
                        {
                            sol.zero();
                            MinimumResidualMethod Solver( stiffness_csr );
//                             PreconditionedConjugateResidualMethod Solver( stiffness_csr, stiffness_invprecon );
                            Solver.print_modulo        = 1+sol.getdimension();
                            Solver.max_iteration_count = 4 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( sol, rhs );
//                             Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2string( end - start ) << std::endl;
                        }


                        cout << "...compute error and residual:" << endl;
            
                        FloatVector error     = interpol_sol  - incmatrix * sol;
                        FloatVector graderror = interpol_grad - diffmatrix * incmatrix * sol;
                        Float errornorm       = std::sqrt( error * ( scalar_massmatrix * error ) );
                        Float graderrornorm   = std::sqrt( graderror * ( vector_massmatrix * graderror ) );
                        Float residualnorm    = ( rhs - stiffness * sol ).norm();
                        
                        cout << "error:     " << errornorm     << endl;
                        cout << "graderror: " << graderrornorm << endl;
                        cout << "residual:  " << residualnorm  << endl;
                        
                        
                        
                        contable << errornorm << graderrornorm << nl;
                        
                        contable.print( std::cout );
                        
                        
                        
                        if( r == 1 ){
                    
                            fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                
                            VTK_MeshWriter_Mesh1D vtk( M, fs );
                            vtk.writePreamble( getbasename(__FILE__) );
                            vtk.writeCoordinateBlock();
                            vtk.writeTopDimensionalCells();
                            
                            vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                            // vtk.writeCellVectorData( interpol_grad, "gradient_interpolation" , 0.1 );
                            
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
