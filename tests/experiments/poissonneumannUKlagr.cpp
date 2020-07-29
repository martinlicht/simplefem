

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
// #include "../../operators/composedoperators.hpp"
#include "../../operators/composedoperators.hpp"
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
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test: Poisson Neumann over " << endl;
        
        cout << std::setprecision(10);

        if(true){

            cout << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = UnitedKingdom(); //
            
            M.check();
            
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            
            
            

            
            // std::function<FloatVector(const std::function<FloatVector(const FloatVector&) ) >scalarfield = 
            
            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector({ 
                        0.0001 * vec[0] * vec[0] * vec[1] * vec[1]
                        });
                };

            std::function<FloatVector(const FloatVector&)> experiment_grad = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector( { 
                            0.0001 * 2. * vec[0] * vec[1] * vec[1],
                            0.0001 * 2. * vec[1] * vec[0] * vec[0], 
                        });
                };

           std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        0.0001 * ( 2. * vec[0] * vec[0] + 2. * vec[1] * vec[1] )
                     });
                };

            

            cout << "Solving Poisson Problem with Neumann boundary conditions over the United Kingdom" << endl;

            int max_l = 2;
            
            ConvergenceTable contable;
            

            for( int l = 0; l <= max_l; l++ ){
                
                cout << "Level: " << l << std::endl;
                cout << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
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
                    auto stiffness = ( incmatrix_t * diffmatrix_t ) * vector_massmatrix * ( diffmatrix * incmatrix );
                    // auto op1_t = incmatrix_t * diffmatrix_t;
                    // auto op2_t = op1_t * vector_massmatrix_fac_t;
                    // auto op1   = diffmatrix * incmatrix;
                    // auto op2   = vector_massmatrix_fac * op1;
                    // auto stiffness = op2_t * op2;

                    // auto opr1 = diffmatrix & incmatrix;
                    // auto opr  = vector_massmatrix_fac & opr1;
                    // auto opl  = opr.getTranspose(); 
                    // auto stiffness = opl & opr;
                    
                    // stiffness.sortentries();
                    // auto stiffness_csr = MatrixCSR( stiffness );
                    
                    auto stiffness_invprecon = DiagonalOperator( stiffness.getdimin(), 1. );
                    // auto stiffness_invprecon = InverseDiagonalPreconditioner( stiffness );
                    // std::cout << "Average value of diagonal preconditioner: " << stiffness_invprecon.getdiagonal().average() << std::endl;

                    {

                        const auto& function_sol  = experiment_sol;
                        const auto& function_grad = experiment_grad;
                        const auto& function_rhs  = experiment_rhs;
                        
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

                        cout << interpol_one * ( scalar_massmatrix * interpol_rhs ) << endl;

                        cout << "...measure interpolation commutativity" << endl;
            
                        cout << "...compute norms of solution and right-hand side:" << endl;
            
                        Float rhs_norm = std::sqrt( interpol_rhs * ( scalar_massmatrix * interpol_rhs ) );

                        cout << "rhs norm:      " << rhs_norm << endl;

                        cout << "...create RHS vector" << endl;

                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( M.count_simplices(0), 0. );
                        
                        cout << "...iterative solver" << endl;
                        
                        sol.zero();
                        FloatVector eins( sol.getdimension(), 1. );
                        
                        sol.random();
                        rhs.zero();
                        
                        if(false)
                        for( int t = 0; t < 4; t++ )
                        {
                            timestamp start = gettimestamp();
                            PreconditionedConjugateResidualMethod CRM( stiffness, stiffness  );
                            CRM.print_modulo = 1;//+sol.getdimension()/1000;
                            //CRM.tolerance = 1e-15;
                            CRM.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << end - start << std::endl;
                            sol = sol - ( interpol_one * ( scalar_massmatrix * incmatrix * sol ) ) / domain_area * eins;
                        }

                        for( int t = 0; t < 4; t++ )
                        {
                            timestamp start = gettimestamp();
                            ConjugateResidualMethod CRM( stiffness );
                            CRM.print_modulo = 1;//+sol.getdimension()/1000;
                            CRM.tolerance = 1e-40;
                            CRM.solve_robust( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << end - start << std::endl;
                            sol = sol - ( interpol_one * ( scalar_massmatrix * incmatrix * sol ) ) / domain_area * eins;
                        }

                        if(false)
                        for( int t = 0; t < 10; t++ )
                        {
                            timestamp start = gettimestamp();
                            ResidualDescentMethod RDM( stiffness );
                            RDM.print_modulo = 1;//+sol.getdimension()/1000;
                            RDM.tolerance = 1e-15;
                            RDM.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << end - start << std::endl;
                            sol = sol - ( interpol_one * ( scalar_massmatrix * incmatrix * sol ) ) / domain_area * eins;
                        }


                        cout << "...compute error and residual:" << endl;
            
                        // Float errornorm     = ( scalar_massmatrix_fac * ( interpol_sol  - incmatrix * sol ) ).norm();
                        // Float graderrornorm = ( vector_massmatrix_fac * ( interpol_grad - diffmatrix * incmatrix * sol ) ).norm();
                        Float residualnorm  = ( rhs - stiffness * sol ).norm();
                        
                        // cout << "error:     " << errornorm     << endl;
                        // cout << "graderror: " << graderrornorm << endl;
                        cout << "residual:  " << residualnorm  << endl;
                        
                        
                        
//                         contable << errornorm << graderrornorm << nl;
//                         
//                         contable.print( std::cout );

                        if( r == 1 ){
                    
                            fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                
                            VTK_MeshWriter_Mesh2D vtk( M, fs );
                            vtk.writePreamble( getbasename(__FILE__) );
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
