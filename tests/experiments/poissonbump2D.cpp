

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
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/minres.hpp"
#include "../../fem/finitediff.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
// #include "../../fem/global.lagrangeincl.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test for Solution of Neumann Problem" << endl;
        
        cout << std::setprecision(10);

        if(true){

            cout << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            M.automatic_dirichlet_flags();
           
            M.check_dirichlet_flags();

            M.getcoordinates().scale(1.1);
            
            cout << "Prepare scalar fields for testing..." << endl;
            


            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector({ 
                        bumpfunction(vec[0]) * bumpfunction(vec[1])
                    });
                };
            
            std::function<FloatVector(const FloatVector&)> experiment_grad = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector( { 
                            bumpfunction_dev(vec[0]) *     bumpfunction(vec[1]),
                            bumpfunction(vec[0])     * bumpfunction_dev(vec[1]), 
                    });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        -
                        bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1])
                        -
                        bumpfunction(vec[0])        * bumpfunction_devdev(vec[1])
                    });
                    

//                     const Float stepsize = 1e-07;
//                     
//                     FloatVector ret(1);
//                     
//                     auto point = vec;
//                     
//                     FloatVector mid    = experiment_sol( point                              );
//                     FloatVector left   = experiment_sol( point - stepsize * unitvector(2,0) );
//                     FloatVector right  = experiment_sol( point + stepsize * unitvector(2,0) );
//                     FloatVector up     = experiment_sol( point - stepsize * unitvector(2,1) );
//                     FloatVector down   = experiment_sol( point + stepsize * unitvector(2,1) );
//                     
//                     ret = - ( up + down + left + right - 4 * mid );
//                     
//                     ret /= ( stepsize * stepsize );
//                     
//                     assert( ret.getdimension() == 1 );
//                     
//                     return ret;
                                    
                };
            
            
            

            

            cout << "Solving Poisson Problem with Dirichlet boundary conditions" << endl;

            int min_l = 1; 
            int max_l = 5;
            
            int min_r = 1;
            int max_r = 1;
            
            ConvergenceTable contable;
            
            contable << "u_error" << "du_error" << nl;
            

            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                cout << "Level: " << l << std::endl;
                cout << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                if( l != 0 )
                for( int r = min_r; r <= max_r; r++ ) 
                {
                    
                    cout << "...assemble scalar mass matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );

                    cout << "...assemble vector mass matrix" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    
                    cout << "...assemble differential matrix and transpose" << endl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    cout << "...assemble inclusion matrix and transpose" << endl;
            
                    SparseMatrix incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );
                    
//                     std::cout << incmatrix.getdimin() <<space<< L_incmatrix.getdimin() <<space<< incmatrix.getdimout() <<space<< L_incmatrix.getdimout() << nl;
//                     assert( incmatrix.getdimin()  == L_incmatrix.getdimin() );
//                     assert( incmatrix.getdimout() == L_incmatrix.getdimout() );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    cout << "...assemble stiffness matrix" << endl;
            
                    // ProductOperator 
//                     auto stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
//                     auto op1 = incmatrix_t * diffmatrix_t;
//                     auto op2 = op1 * vector_massmatrix;
//                     auto op3 = op2 * diffmatrix;
//                     auto stiffness = op3 * incmatrix;
//                     auto& stiffness_csr = stiffness;

                    auto opr  = diffmatrix & incmatrix;
                    auto opl  = opr.getTranspose(); 
                    auto stiffness = opl & ( vector_massmatrix & opr );
                    
                    stiffness.sortentries();
                    auto stiffness_csr = MatrixCSR( stiffness );
                    
                    auto stiffness_invprecon = DiagonalOperator( stiffness.getdimin(), 1. );
//                     auto stiffness_invprecon = InverseDiagonalPreconditioner( stiffness );
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
            
                        auto commutatorerror_aux = interpol_grad - diffmatrix * interpol_sol;
                        Float commutatorerror = commutatorerror_aux * ( vector_massmatrix * commutatorerror_aux );
                        cout << "commutator error: " << commutatorerror << endl;
                        
                        cout << "...compute norms of solution and right-hand side:" << endl;
            
                        Float sol_norm = interpol_sol * ( scalar_massmatrix * interpol_sol );
                        Float rhs_norm = interpol_rhs * ( scalar_massmatrix * interpol_rhs );
                        
                        cout << "solution norm: " << sol_norm << endl;
                        cout << "rhs norm:      " << rhs_norm << endl;

                        cout << "...create RHS vector" << endl;

                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( incmatrix.getdimin(), 0. );
                        
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
                            std::cout << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        }

                        cout << "...compute error and residual:" << endl;
            
                        
                        auto errornorm_aux = interpol_sol  - incmatrix * sol;
                        auto graderror_aux = interpol_grad - diffmatrix * incmatrix * sol;
                        
                        Float errornorm     = sqrt( errornorm_aux * ( scalar_massmatrix * errornorm_aux ) );
                        Float graderrornorm = sqrt( graderror_aux * ( vector_massmatrix * graderror_aux ) );
                        Float residualnorm  = ( rhs - stiffness * sol ).norm();
                        
                        cout << "error:     " << errornorm     << endl;
                        cout << "graderror: " << graderrornorm << endl;
                        cout << "residual:  " << residualnorm  << endl;
                        
                        
                        
                        contable << errornorm << graderrornorm << nl;
                        
                        contable.print( std::cout );


                        if( r == 1 ){
                    
                            
                            auto outputdata1 = sol;
                            auto outputdata2 = sol;
                            
                            
                            for( int c = 0; c < M.count_simplices(0); c++ ) { 
                                auto x = M.getcoordinates().getdata(c,0);
                                auto y = M.getcoordinates().getdata(c,1);
                                auto value = experiment_rhs( { x, y } )[0];
                                outputdata2[c] = value;
                            }
                            
                            
                            fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                
                            VTKWriter vtk( M, fs, getbasename(__FILE__) );
                            vtk.writeCoordinateBlock( outputdata1 );
                            vtk.writeTopDimensionalCells();
                            
                            vtk.writeVertexScalarData( -outputdata1, "iterativesolution_scalar_data" , 1.0 );
                            vtk.writeVertexScalarData(  outputdata2, "reference_scalar_data" , 1.0 );
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
