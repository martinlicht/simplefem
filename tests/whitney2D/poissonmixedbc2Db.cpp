

/**/

#include <ostream>
#include <fstream>
// #include <iomanip>

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
#include "../../solver/sparsesolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test for Solution of Neumann Problem" << endl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            M.set_flag( 1, 0, SimplexFlagDirichlet );
            M.set_flag( 1, 1, SimplexFlagDirichlet );
            M.set_flag( 0, 0, SimplexFlagDirichlet );
            M.set_flag( 0, 1, SimplexFlagDirichlet );
            M.set_flag( 0, 2, SimplexFlagDirichlet );

            
            LOG << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            
            
            


            
            // std::function<FloatVector(const std::function<FloatVector(const FloatVector&) ) >scalarfield = 
            
            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    Float x = vec[0]; Float y = vec[1];
                    Float k = 2*Constants::pi;
                    return FloatVector({ square( sin(k*x) * sin(k*y) ) });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_grad = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    Float x = vec[0]; Float y = vec[1];
                    Float k = 2*Constants::pi;
                    return FloatVector( { 
                        k * sin( 2*k*x ) * sin(k*y) * sin(k*y),
                        k * sin( 2*k*y ) * sin(k*x) * sin(k*x),
                    });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    Float x  =  vec[0]; Float y  =  vec[1];
                    Float k = 2*Constants::pi;
                    return FloatVector({ 
                        - k*k * ( cos(2*k*x) - cos( 2*k*(x - y) ) + cos(2*k*y) - cos( 2*k*(x + y) ) ) 
                     });
                };
            

            

            

            LOG << "Solving Poisson Problem with Neumann boundary conditions" << endl;

            const int min_l = 0; 
            const int max_l = 8;
            
            const int min_r = 3;
            const int max_r = 3;
            
            ConvergenceTable contable("Mass error");
            
            contable << "u_error" << "du_error" << "residual" << "time" << nl;
            

            assert( 0 <= min_l and min_l <= max_l );
            assert( 0 <= min_r and min_r <= max_r );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << "/" << max_l << std::endl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                if( l != 0 )
                for( int r = min_r; r <= max_r; r++ ) 
                {
                    
                    LOG << "...assemble scalar mass matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );

                    LOG << "...assemble vector mass matrix" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    
                    LOG << "...assemble differential matrix and transpose" << endl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    LOG << "...assemble inclusion matrix and transpose" << endl;
            
                    SparseMatrix incmatrix = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r );
                    
//                     LOG << incmatrix.getdimin() <<space<< L_incmatrix.getdimin() <<space<< incmatrix.getdimout() <<space<< L_incmatrix.getdimout() << nl;
//                     assert( incmatrix.getdimin()  == L_incmatrix.getdimin() );
//                     assert( incmatrix.getdimout() == L_incmatrix.getdimout() );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    LOG << "...assemble stiffness matrix" << endl;
            
                    // ProductOperator 
//                     auto stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
//                     auto op1 = incmatrix_t * diffmatrix_t;
//                     auto op2 = op1 * vector_massmatrix;
//                     auto op3 = op2 * diffmatrix;
//                     auto stiffness = op3 * incmatrix;
//                     auto& stiffness_csr = stiffness;

                    auto opr = diffmatrix & incmatrix;
                    auto opl = opr.getTranspose(); 
                    auto stiffness = opl & ( vector_massmatrix & opr );
                    
                    stiffness.sortentries();
                    auto stiffness_csr = MatrixCSR( stiffness );
                    
                    auto stiffness_invprecon = DiagonalOperator( stiffness.getdimin(), 1. );
//                     auto stiffness_invprecon = InverseDiagonalPreconditioner( stiffness );
                    LOG << "Average value of diagonal preconditioner: " << stiffness_invprecon.getdiagonal().average() << std::endl;

                    {

                        const auto& function_sol  = experiment_sol;
                        const auto& function_grad = experiment_grad;
                        const auto& function_rhs  = experiment_rhs;
                        
                        LOG << "...interpolate explicit solution and rhs" << endl;
            
                        FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, r,   function_sol  );
                        FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r-1, function_grad );
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        
                        LOG << "...compute norms of solution and right-hand side:" << endl;
            
                        Float sol_norm = interpol_sol * ( scalar_massmatrix * interpol_sol );
                        Float rhs_norm = interpol_rhs * ( scalar_massmatrix * interpol_rhs );
                        
                        LOG << "solution norm: " << sol_norm << endl;
                        LOG << "rhs norm:      " << rhs_norm << endl;

                        LOG << "...create RHS vector" << endl;

                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( incmatrix.getdimin(), 0. );
                        
                        LOG << "...iterative solver" << endl;
                        
                        timestamp start = gettimestamp();

                        {
                            sol.zero();
                            FloatVector residual( rhs );
                            
                            ConjugateGradientSolverCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness_csr.getA(), stiffness_csr.getC(), stiffness_csr.getV(),
                                residual.raw(),
                                1e-16,
                                1
                            );

                        }

                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                        LOG << "...compute error and residual:" << endl;
            
                        
                        auto errornorm_aux = interpol_sol  - incmatrix * sol;
                        auto graderror_aux = interpol_grad - diffmatrix * incmatrix * sol;
                        
                        Float errornorm     = sqrt( errornorm_aux * ( scalar_massmatrix * errornorm_aux ) );
                        Float graderrornorm = sqrt( graderror_aux * ( vector_massmatrix * graderror_aux ) );
                        Float residualnorm  = ( rhs - stiffness * sol ).norm();
                        
                        LOG << "error:     " << errornorm    << endl;
                        LOG << "graderror: " << graderrornorm << endl;
                        LOG << "residual:  " << residualnorm << endl;
                        LOG << "time:      " << Float( end - start ) << endl;
                        
                        contable << errornorm;
                        contable << graderrornorm;
                        contable << residualnorm;
                        contable << Float( end - start );
                        contable << nl;
                        
                        contable.lg();


                        if( r == 1 ){
                        
                            fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                        
                            VTKWriter vtk( M, fs, getbasename(__FILE__) );
                            vtk.writeCoordinateBlock( 0.3 * sol );
                            vtk.writeTopDimensionalCells();
                            
                            vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                            // vtk.writeCellVectorData( interpol_grad, "gradient_interpolation" , 0.1 );
                            
                            fs.close();
                    
                        }


                    }
                    
                }

                if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
                
                

            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
