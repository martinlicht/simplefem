

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
// #include "../../operators/composed.hpp"
#include "../../dense/densematrix.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../vtk/vtkwriter.mesh1D.hpp"
#include "../../solver/crm.hpp"
#include "../../solver/pcrm.hpp"
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

            cout << "Case 1D" << endl;
            
            cout << "Initial mesh..." << endl;
            
            MeshSimplicial1D M = StandardSquare1D();
            
            M.check();
            
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 1 );
                        return FloatVector({ 1. });
                    };
            
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_rhs;
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_grad;
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_sol;


            
            // std::function<FloatVector(const FloatVector&) scalarfield = 
            
            Float xfeq = 1.;
            
            experiments_sol.push_back( 
                [xfeq](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 1 );
                    return FloatVector({ std::cos( xfeq * Constants::twopi * vec[0] ) });
                }
            );

            experiments_grad.push_back( 
                [xfeq](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 1 );
                    return FloatVector( { 
                            - xfeq * Constants::twopi * std::sin( xfeq * Constants::twopi * vec[0] )
                        });
                }
            );

            experiments_rhs.push_back( 
                [xfeq](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 1 );
                    return FloatVector({ 
                        xfeq*xfeq * Constants::fourpisquare * std::cos( xfeq * Constants::twopi * vec[0] )
                     });
                }
            );

            

            ConvergenceTable contable;
            

            assert( experiments_sol.size() == experiments_rhs.size() && experiments_sol.size() == experiments_grad.size() );

            cout << "Solving Poisson Problem with Dirichlet boundary conditions" << endl;

            int max_l = 15;
            
            for( int l = 0; l <= max_l; l++ ){
                
                cout << "Level: " << l << std::endl;
                cout << "# E/V: " << M.count_edges() << "/" << M.count_vertices() << nl;
                
                const int r = 1;
                
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

                    auto opr1 = diffmatrix & incmatrix;
                    auto opr  = vector_massmatrix_fac & opr1;
                    auto opl  = opr.getTranspose(); 
                    auto stiffness = opl & opr;
                    
                    stiffness.sortentries();
                    auto stiffness_csr = MatrixCSR( stiffness );
                    
                    auto stiffness_invprecon = DiagonalOperator( stiffness.getdimin(), 1. );
                    //auto stiffness_invprecon = InverseDiagonalPreconditioner( stiffness );
                    std::cout << "Average value of diagonal preconditioner: " << stiffness_invprecon.getdiagonal().average() << std::endl;

                    for( int i = 0; i < experiments_sol.size(); i++){

                        const auto& function_sol = experiments_sol[i];
                        const auto& function_grad= experiments_grad[i];
                        const auto& function_rhs = experiments_rhs[i];
                        
                        cout << "...interpolate explicit solution and rhs" << endl;
            
                        FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, r,   function_sol  );
                        FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r-1, function_grad );
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        
                        FloatVector interpol_one  = Interpolation( M, M.getinnerdimension(), 0, r, constant_one );
                        
                        cout << "...measure kernel component: " << std::flush;
            
                        Float average_sol = interpol_one * ( scalar_massmatrix * interpol_sol );
                        Float average_rhs = interpol_one * ( scalar_massmatrix * interpol_rhs );
                        
                        cout << average_sol << space << average_rhs << endl;

                        cout << "...measure interpolation commutativity" << endl;
            
                        Float commutatorerror = ( vector_massmatrix_fac * ( interpol_grad - diffmatrix * interpol_sol ) ).norm();
                        cout << "commutator error: " << commutatorerror << endl;
                        
                        cout << "...compute norms of solution and right-hand side:" << endl;
            
                        Float sol_norm = ( scalar_massmatrix_fac * interpol_sol ).norm();
                        Float rhs_norm = ( scalar_massmatrix_fac * interpol_rhs ).norm();
                        
                        cout << "solution norm: " << sol_norm << endl;
                        cout << "rhs norm:      " << rhs_norm << endl;

                        cout << "...create RHS vector" << endl;

                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( M.count_simplices(0), 0. );
                        
                        cout << "...iterative solver" << endl;
                        
                        {
                            sol.zero();
                            timestamp start = gettimestamp();
                            ConjugateResidualMethod CRM( stiffness_csr );
                            CRM.print_modulo = 1+sol.getdimension();
                            CRM.tolerance = 1e-50;
                            CRM.solve_robust( sol, rhs );
                            CRM.solve_robust( sol, rhs );
                            CRM.solve_robust( sol, rhs );
                            CRM.solve_robust( sol, rhs );
                            CRM.solve_robust( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << end - start << std::endl;
                        }
                        
                        if(false){
                            sol.zero();
                            timestamp start = gettimestamp();
                            PreconditionedConjugateResidualMethod PCRM( stiffness_csr, stiffness_invprecon );
                            PCRM.print_modulo = 1+sol.getdimension();
                            PCRM.tolerance = 1e-10;
                            PCRM.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << end - start << std::endl;
                        }

                        cout << "...compute error and residual:" << endl;
            
                        Float average_numsol = interpol_one * ( scalar_massmatrix * incmatrix * sol );
                        
                        Float errornorm     = ( scalar_massmatrix_fac * ( interpol_sol  - incmatrix * sol - average_numsol * interpol_one ) ).norm();
                        Float graderrornorm = ( vector_massmatrix_fac * ( interpol_grad - diffmatrix * incmatrix * sol ) ).norm();
                        Float residualnorm  = ( rhs - stiffness * sol ).norm();
                        
                        // FloatVector gradfoo = diffmatrix * ( interpol_sol - incmatrix * sol );
                        // Float graderrornorm = gradfoo.scalarproductwith( vector_massmatrix * gradfoo );
                        // Float errornorm1 = interpol_sol * ( scalar_massmatrix * interpol_sol );
                        // Float errornorm2 = power( ( scalar_massmatrix_fac * interpol_sol ).norm(), 2. );

                        cout << "error:     " << errornorm     << endl;
                        cout << "graderror: " << graderrornorm << endl;
                        cout << "residual:  " << residualnorm  << endl;
                        
                        
                        
                        contable << errornorm << graderrornorm << average_numsol << nl;
                        
                        contable.print( std::cout );
                        
                        
                        
//                         {
//                     
//                             fstream fs( adaptfilename("./poissondirichlet.vtk"), std::fstream::out );
//                 
//                             VTK_MeshWriter_Mesh1D vtk( M, fs );
//                             vtk.writePreamble( "Poisson-Dirichlet problem" );
//                             vtk.writeCoordinateBlock( 0.3 * sol );
//                             vtk.writeTopDimensionalCells();
//                             
//                             vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
//                             // vtk.writeCellVectorData( interpol_grad, "gradient_interpolation" , 0.1 );
//                             
//                             fs.close();
//                     
//                         }


                    }
                    
                }

                cout << "Refinement..." << endl;
            
                if( l != max_l ) M.uniformrefinement();
                
                

            } 
        
        }
        
        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
