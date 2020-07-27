

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
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.mesh2D.hpp"
#include "../../solver/crm.hpp"
#include "../../solver/pcrm.hpp"
#include "../../solver/minres.hpp"
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
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_grad;
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

            experiments_grad.push_back( 
                [xfeq,yfeq](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector( { 
                            -xfeq * Constants::pi * std::sin( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ),
                            -yfeq * Constants::pi * std::cos( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] ), 
                        });
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

            int max_l = 8;
            int max_r = 5;
            
            for( int l = 0; l <= max_l; l++ ){
                
                cout << "Level: " << l << std::endl;
                cout << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                for( int r = 5; r <= max_r; r++ ) 
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
                    auto op1 = incmatrix_t * diffmatrix_t;
                    auto op2 = op1 * vector_massmatrix;
                    auto op3 = op2 * diffmatrix;
                    auto stiffness = op3 * incmatrix;
                    auto& stiffness_csr = stiffness;

//                     auto opr = diffmatrix & incmatrix;
//                     auto opl = opr.getTranspose(); 
//                     auto stiffness = opl & ( vector_massmatrix & opr );
                    
//                     stiffness.sortentries();
//                     auto stiffness_csr = MatrixCSR( stiffness );
                    
                    auto stiffness_invprecon = DiagonalOperator( stiffness.getdimin(), 1. );
//                     auto stiffness_invprecon = InverseDiagonalPreconditioner( stiffness );
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
            
                        auto commutatorerror_aux = interpol_grad - diffmatrix * interpol_sol;
                        Float commutatorerror = commutatorerror_aux * ( vector_massmatrix * commutatorerror_aux );
                        cout << "commutator error: " << commutatorerror << endl;
                        
                        cout << "...compute norms of solution and right-hand side:" << endl;
            
                        Float sol_norm = interpol_sol * ( scalar_massmatrix * interpol_sol );
                        Float rhs_norm = interpol_sol * ( scalar_massmatrix * interpol_rhs );
                        
                        cout << "solution norm: " << sol_norm << endl;
                        cout << "rhs norm:      " << rhs_norm << endl;

                        cout << "...create RHS vector" << endl;

                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( incmatrix.getdimin(), 0. );
                        
                        cout << "...iterative solver" << endl;
                        
                        {
                            sol.zero();
                            timestamp start = gettimestamp();
                            MinimumResidualMethod MINRES( stiffness_csr );
                            MINRES.print_modulo = 1+sol.getdimension();
                            MINRES.tolerance = 1e-20;
//                             MINRES.solve_robust( sol, rhs );
                            MINRES.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t " << end - start << std::endl;
                        }
                        
                        if(false){
                            sol.zero();
                            timestamp start = gettimestamp();
                            PreconditionedConjugateResidualMethod PCRM( stiffness_csr, stiffness_invprecon );
                            PCRM.print_modulo = 1+sol.getdimension();
                            PCRM.tolerance = 1e-10;
                            PCRM.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t " << end - start << std::endl;
                        }

                        cout << "...compute error and residual:" << endl;
            
                        
                        auto errornorm_aux = interpol_sol  - incmatrix * sol;
                        auto graderror_aux = interpol_grad - diffmatrix * incmatrix * sol;
                        
                        Float errornorm     = sqrt( errornorm_aux * ( scalar_massmatrix * errornorm_aux ) );
                        Float graderrornorm = sqrt( graderror_aux * ( vector_massmatrix * graderror_aux ) );
                        Float residualnorm  = ( rhs - stiffness * sol ).norm();
                        
                        // FloatVector gradfoo = diffmatrix * ( interpol_sol - incmatrix * sol );
                        // Float graderrornorm = gradfoo.scalarproductwith( vector_massmatrix * gradfoo );
                        // Float errornorm1 = interpol_sol * ( scalar_massmatrix * interpol_sol );
                        // Float errornorm2 = power( ( scalar_massmatrix_fac * interpol_sol ).norm(), 2. );

                        cout << "error:     " << errornorm     << endl;
                        cout << "graderror: " << graderrornorm << endl;
                        cout << "residual:  " << residualnorm  << endl;


                        if( r == 1 ){
                    
                            fstream fs( adaptfilename("./poissonneumann.vtk"), std::fstream::out );
                
                            VTK_MeshWriter_Mesh2D vtk( M, fs );
                            vtk.writePreamble( "Poisson-Neumann problem" );
                            vtk.writeCoordinateBlock( 0.3 * sol );
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
