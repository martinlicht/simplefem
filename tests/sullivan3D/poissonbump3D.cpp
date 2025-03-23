

/**/

#include <cmath>

#include <fstream>
#include <functional>


#include "../../base/include.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../utility/math.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.inclsullivan.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: 3D Poisson Problem" << nl;
        
    if(true){

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial3D M = StandardCube3D();
        
        M.check();
        
        M.automatic_dirichlet_flags();
        
        M.check_dirichlet_flags();

        M.getCoordinates().scale(1.1);
        
        LOG << "Prepare scalar fields for testing..." << nl;
        


        std::function<FloatVector(const FloatVector&)> experiment_sol = 
            [=](const FloatVector& vec) -> FloatVector {
                assert( vec.getdimension() == 3 );
                // return FloatVector({ 1. });
                return FloatVector({ 
                    bumpfunction(vec[0]) * bumpfunction(vec[1]) * bumpfunction(vec[2])
                });
            };
        
        std::function<FloatVector(const FloatVector&)> experiment_grad = 
            [=](const FloatVector& vec) -> FloatVector {
                assert( vec.getdimension() == 3 );
                return FloatVector( { 
                        bumpfunction_dev(vec[0]) *     bumpfunction(vec[1]) *     bumpfunction(vec[2]),
                            bumpfunction(vec[0]) * bumpfunction_dev(vec[1]) *     bumpfunction(vec[2]), 
                            bumpfunction(vec[0]) *     bumpfunction(vec[1]) * bumpfunction_dev(vec[2]), 
                });
            };
        

        std::function<FloatVector(const FloatVector&)> experiment_rhs = 
            [=](const FloatVector& vec) -> FloatVector {
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
        const int max_l = 3;
        
        const int min_r = 1;
        const int max_r = 1;
        
        ConvergenceTable contable("Mass error");
        
        contable << "u_error" << "du_error" << "residual" << "time" << nl;
        

        assert( 0 <= min_l and min_l <= max_l );
        assert( 0 <= min_r and min_r <= max_r );
        
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ ){
            
            LOG << "Level: " << min_l << " <= " << l << " <= " << max_l << nl;
            LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            if( l != 0 )
            for( int r = min_r; r <= max_r; r++ ) 
            {
                
                LOG << "Level: "             << min_l << " <= " << l << " <= " << max_l << nl;
                LOG << "Polynomial degree: " << min_r << " <= " << r << " <= " << max_r << nl;
                    
                LOG << "... assemble scalar mass matrices" << nl;
        
                SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );

                LOG << "... assemble vector mass matrix" << nl;
        
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                
                LOG << "... assemble inclusion matrix and transpose" << nl;
        
                SparseMatrix incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );
                
                SparseMatrix incmatrix_t = incmatrix.getTranspose();

                LOG << "... assemble differential matrix and transpose" << nl;

                SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                LOG << "... assemble stiffness matrix" << nl;

                // auto opr  = diffmatrix & incmatrix;
                // auto opl  = opr.getTranspose(); 
                // auto stiffness = opl & ( vector_massmatrix & opr );
                // stiffness.sortentries();
                // auto stiffness_csr = MatrixCSR( stiffness );
                
                auto stiffness_csr = Conjugation( MatrixCSR(vector_massmatrix), MatrixCSR(diffmatrix) & MatrixCSR(incmatrix) );

                // LOG << "... compose preconditioner data" << nl;

                {

                    LOG << "... interpolate explicit solution and rhs" << nl;
        
                    const auto& function_sol  = experiment_sol;
                    const auto& function_grad = experiment_grad;
                    const auto& function_rhs  = experiment_rhs;
                    
                    FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, r,   function_sol  );
                    FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r-1, function_grad );
                    FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                    
                    LOG << "... compute norms of solution and right-hand side" << nl;
        
                    Float sol_norm = interpol_sol * ( scalar_massmatrix * interpol_sol );
                    Float rhs_norm = interpol_rhs * ( scalar_massmatrix * interpol_rhs );
                    
                    LOG << "solution norm: " << sol_norm << nl;
                    LOG << "rhs norm:      " << rhs_norm << nl;

                    FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                    FloatVector sol( incmatrix.getdimin(), 0. );
                    
                    LOG << "... iterative solver" << nl;
                    
                    timestamp start = timestampnow();
                    
                    {
                        sol.zero();
                        ConjugateGradientMethod solver( stiffness_csr );
                        solver.solve( sol, rhs );
                    }

                    timestamp end = timestampnow();
                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                    LOG << "... compute error and residual" << nl;
        
                    
                    auto computed_sol  = incmatrix * sol;
                    auto computed_grad = diffmatrix * incmatrix * sol;
                    
                    auto errornorm_aux = interpol_sol  - computed_sol;
                    auto graderror_aux = interpol_grad - computed_grad;
                    
                    Float errornorm     = std::sqrt( errornorm_aux * ( scalar_massmatrix * errornorm_aux ) );
                    Float graderrornorm = std::sqrt( graderror_aux * ( vector_massmatrix * graderror_aux ) );
                    Float residualnorm  = ( rhs - stiffness_csr * sol ).norm();
                    
                    LOG << "error:     " << errornorm    << nl;
                    LOG << "graderror: " << graderrornorm << nl;
                    LOG << "residual:  " << residualnorm << nl;
                    
                    
                    
                    contable << errornorm;
                    contable << graderrornorm;
                    contable << residualnorm;
                    contable << Float( end - start );
                    contable << nl;
                    
                    contable.lg();


                    {
                        /* Print the iterative solution and also the data */

                        std::fstream fs( get_available_filename(get_basename(__FILE__)), std::fstream::out );
                        VTKWriter vtk( M, fs, get_basename(__FILE__) );
                        
                        vtk.write_vertex_scalar_data( [&]( const FloatVector& vec) -> Float{ return function_sol(vec)[0]; }, "interpolated_vertex_sol" );
                        vtk.write_vertex_scalar_data( [&]( const FloatVector& vec) -> Float{ return function_rhs(vec)[0]; }, "interpolated_vertex_rhs" );

                        if( r == 1) {
                            vtk.write_vertex_scalar_data( sol, "iterativesolution_scalar_data" , 1.0 );
                        }
                        
                        vtk.write_cell_scalar_data( [&]( const FloatVector& vec) -> Float{ return function_sol(vec)[0]; }, "interpolated_sol" );
                        vtk.write_cell_scalar_data( [&]( const FloatVector& vec) -> Float{ return function_rhs(vec)[0]; }, "interpolated_rhs" );

                        {
                            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r );
                            const auto printable_sol = interpol_matrix * incmatrix * sol; 
                            vtk.write_cell_scalar_data( printable_sol, "iterativesolution_scalar_data_cellwise" , 1.0 );
                        }
                        

                        {
                            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r-1 );
                            const auto printable_grad = interpol_matrix * diffmatrix * incmatrix * sol;
                        
                            vtk.write_cell_vector_data_barycentricgradients( printable_grad,  "gradient_calculation" );
                            vtk.write_cell_vector_data( function_grad,  "gradient_interpolation" );
                        }
                        
                        {
                            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r-1 );
                            const auto printable_grad = interpol_matrix * computed_grad; 
                            vtk.write_cell_vector_data_barycentricgradients( printable_grad, "gradient_interpolation" , 1.0 );
                        }
                        
                        fs.close();
                    }


                }
                
            }

            if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
    
        } 
    
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
