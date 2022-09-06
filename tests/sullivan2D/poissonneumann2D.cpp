

/**/

#include <fstream>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test for Solution of Neumann Problem" << nl;
        
        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            LOG << "Prepare scalar fields for testing..." << nl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            
            
            


            
            const Float xfeq = 1.;
            const Float yfeq = 1.;
            

            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector({ std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ) });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_grad = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector( { 
                            -xfeq * Constants::pi * std::sin( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ),
                            -yfeq * Constants::pi * std::cos( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] ), 
                        });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        xfeq*xfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
                        +
                        yfeq*yfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
                     });
                };
            

            

            

            LOG << "Solving Poisson Problem with Neumann boundary conditions" << nl;

            const int min_l = 0; 
            const int max_l = 8;
            
            const int min_r = 1;
            const int max_r = 1;
            
            const int r_plus_scalar = 0;
            const int r_plus_vector = 0;
            
            ConvergenceTable contable("Mass error");
            
            contable << "u_error" << "du_error" << "residual" << "time" << nl;
            

            assert( 0 <= min_l and min_l <= max_l );
            assert( 0 <= min_r and min_r <= max_r );
            assert( 0 <= r_plus_scalar and 0 <= r_plus_vector);
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << "/" << max_l << nl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                for( int r = min_r; r <= max_r; r++ ) 
                {
                    
                    LOG << "...assemble scalar mass matrices" << nl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r   + r_plus_scalar );

                    LOG << "...assemble vector mass matrix" << nl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 + r_plus_vector );
                    
                    LOG << "...assemble differential matrix and transpose" << nl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    LOG << "...assemble inclusion matrix and transpose" << nl;
            
                    SparseMatrix incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );
                    
                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    LOG << "...assemble elevation matrices and transposes" << nl;
            
                    SparseMatrix scalar_elevmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r  , r_plus_scalar );
                    SparseMatrix vector_elevmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r-1, r_plus_vector );
                    
                    SparseMatrix scalar_elevmatrix_t = scalar_elevmatrix.getTranspose();
                    SparseMatrix vector_elevmatrix_t = vector_elevmatrix.getTranspose();

                    LOG << "...assemble stiffness matrix" << nl;
            
                    // ProductOperator 
//                     auto stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
                    auto op1 = incmatrix_t * diffmatrix_t * vector_elevmatrix_t;
                    auto op2 = op1 * vector_massmatrix;
                    auto op3 = op2 * vector_elevmatrix * diffmatrix;
                    auto stiffness = op3 * incmatrix;
                    auto& stiffness_csr = stiffness;

//                     auto opr  = diffmatrix & incmatrix;
//                     auto opl  = opr.getTranspose(); 
//                     auto stiffness = opl & ( vector_massmatrix & opr );
                    
//                     stiffness.sortentries();
//                     auto stiffness_csr = MatrixCSR( stiffness );
                    
                    auto stiffness_invprecon = DiagonalOperator( stiffness.getdimin(), 1. );
//                     auto stiffness_invprecon = InverseDiagonalPreconditioner( stiffness );
                    LOG << "Average value of diagonal preconditioner: " << stiffness_invprecon.getdiagonal().average() << nl;

                    {

                        const auto& function_sol  = experiment_sol;
                        const auto& function_grad = experiment_grad;
                        const auto& function_rhs  = experiment_rhs;
                        
                        LOG << "...interpolate explicit solution and rhs" << nl;
            
                        FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, r   + r_plus_scalar, function_sol  );
                        FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r-1 + r_plus_vector, function_grad );
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r   + r_plus_scalar, function_rhs  );
                        
                        FloatVector interpol_one  = Interpolation( M, M.getinnerdimension(), 0, r, constant_one );
                        
                        LOG << "...measure kernel component: ";
            
                        Float average_sol = interpol_one * ( scalar_massmatrix * interpol_sol );
                        Float average_rhs = interpol_one * ( scalar_massmatrix * interpol_rhs );
                        
                        LOG << average_sol << space << average_rhs << nl;

                        LOG << "...compute norms of solution and right-hand side:" << nl;
            
                        Float sol_norm = interpol_sol * ( scalar_massmatrix * interpol_sol );
                        Float rhs_norm = interpol_rhs * ( scalar_massmatrix * interpol_rhs );
                        
                        LOG << "solution norm: " << sol_norm << nl;
                        LOG << "rhs norm:      " << rhs_norm << nl;

                        LOG << "...create RHS vector" << nl;

                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( incmatrix.getdimin(), 0. );
                        
                        LOG << "...iterative solver" << nl;
                        
                        timestamp start = gettimestamp();

                        {
                            sol.zero();
                            ConjugateResidualMethod Solver( stiffness_csr );
                            Solver.solve( sol, rhs );
                        }

                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                        LOG << "...compute error and residual:" << nl;
            
                        
                        auto computed_sol  = incmatrix * sol;
                        auto computed_grad = diffmatrix * incmatrix * sol;

                        auto polluted_error = interpol_sol  - scalar_elevmatrix * incmatrix * sol;
                        Float average_error = ( interpol_one * ( scalar_massmatrix * polluted_error ) ) / ( interpol_one * ( scalar_massmatrix * interpol_one ) );
                        auto          error = polluted_error - average_error * interpol_one;
                        auto      graderror = interpol_grad - vector_elevmatrix * diffmatrix * incmatrix * sol;
                        
                        Float errornorm       = std::sqrt( error * ( scalar_massmatrix * error ) );
                        Float graderrornorm   = std::sqrt( graderror * ( vector_massmatrix * graderror ) );
                        Float residualnorm    = ( rhs - stiffness * sol ).norm();
                        
                        LOG << "error:     " << errornorm    << nl;
                        LOG << "graderror: " << graderrornorm << nl;
                        LOG << "residual:  " << residualnorm << nl;
                        LOG << "time:      " << Float( end - start ) << nl;
                        
                        contable << errornorm;
                        contable << graderrornorm;
                        contable << residualnorm;
                        contable << Float( end - start );
                        contable << nl;
                        
                        contable.lg();


                        if( r == 1 ){
                            fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                            VTKWriter vtk( M, fs, getbasename(__FILE__) );
                            vtk.writeCoordinateBlock();
                            vtk.writeTopDimensionalCells();
                            vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                            vtk.writeCellVectorData( computed_grad, "gradient_interpolation" , 0.1 );
                            fs.close();
                        }


                    }
                    
                }

                if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
                
                

            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}
