

/**/

#include <fstream>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/minres.hpp"
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
            
            MeshSimplicial3D M = StandardCube3D();
            
            M.check();
            
            M.automatic_dirichlet_flags();
           
            M.check_dirichlet_flags();

            M.getcoordinates().scale(1.1);
            
            LOG << "Prepare scalar fields for testing..." << nl;
            


            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    // return FloatVector({ 1. });
                    return FloatVector({ 
                        bumpfunction(vec[0]) * bumpfunction(vec[1]) * bumpfunction(vec[2])
                    });
                };
            
            std::function<FloatVector(const FloatVector&)> experiment_grad = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    return FloatVector( { 
                            bumpfunction_dev(vec[0]) *     bumpfunction(vec[1]) *     bumpfunction(vec[2]),
                                bumpfunction(vec[0]) * bumpfunction_dev(vec[1]) *     bumpfunction(vec[2]), 
                                bumpfunction(vec[0]) *     bumpfunction(vec[1]) * bumpfunction_dev(vec[2]), 
                    });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    return FloatVector({ 
                        -
                        bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1]) *        bumpfunction(vec[2])
                        -
                        bumpfunction(vec[0])        * bumpfunction_devdev(vec[1]) *        bumpfunction(vec[2])
                        -
                        bumpfunction(vec[0])        *        bumpfunction(vec[1]) * bumpfunction_devdev(vec[2])
                    });
                
                };
            
            
            

            

            const int min_l = 0; 
            const int max_l = 4;
            
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
                LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                if( l != 0 )
                for( int r = min_r; r <= max_r; r++ ) 
                {
                    
                    LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                    
                    LOG << "integration with: " << r_plus_scalar << ", " << r_plus_vector << nl;
                    
                    LOG << "...assemble scalar mass matrices" << nl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );

                    LOG << "...assemble vector mass matrix" << nl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    
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

                    auto opr  = diffmatrix & incmatrix;
                    auto opl  = opr.getTranspose(); 
                    auto stiffness = opl & ( vector_massmatrix & opr );
                    
                    stiffness.sortentries();
                    auto stiffness_csr = MatrixCSR( stiffness );
                    
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
                        
                        // LOG << "...compute norms of solution and right-hand side:" << nl;
            
                        // Float sol_norm = interpol_sol * ( scalar_massmatrix * interpol_sol );
                        // Float rhs_norm = interpol_rhs * ( scalar_massmatrix * interpol_rhs );
                        
                        // LOG << "solution norm: " << sol_norm << nl;
                        // LOG << "rhs norm:      " << rhs_norm << nl;

                        LOG << "...create RHS vector" << nl;

                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( incmatrix.getdimin(), 0. );
                        
                        LOG << "...iterative solver" << nl;
                        
                        timestamp start = gettimestamp();
                        
                        {
                            sol.zero();
                            ConjugateGradientMethod Solver( stiffness_csr );
                            Solver.solve( sol, rhs );
                        }

                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                        LOG << "...compute error and residual:" << nl;
            
                        
                        auto computed_sol  = incmatrix * sol;
                        auto computed_grad = diffmatrix * incmatrix * sol;
                        
                        auto errornorm_aux = interpol_sol  - scalar_elevmatrix * computed_sol;
                        auto graderror_aux = interpol_grad - vector_elevmatrix * computed_grad;
                        
                        Float errornorm     = sqrt( errornorm_aux * ( scalar_massmatrix * errornorm_aux ) );
                        Float graderrornorm = sqrt( graderror_aux * ( vector_massmatrix * graderror_aux ) );
                        Float residualnorm  = ( rhs - stiffness * sol ).norm();
                        
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
