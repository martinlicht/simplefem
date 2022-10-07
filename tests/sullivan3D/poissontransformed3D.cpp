

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
#include "../../solver/sparsesolver.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.coefficientmassmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
    LOG << "Unit Test for Solution of Neumann Problem" << nl;
    
    if(true){

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial3D M = FicheraCorner3D();
        
        M.check();
        
        M.automatic_dirichlet_flags();
        
        for( int d = 0; d <= 2; d++ )
        for( int s = 0; s < M.count_simplices(d); s++ )
            if( M.get_flag(d,s) == SimplexFlagNull and M.get_midpoint(d,s)[0] >= 0. )
                M.set_flag( d, s, SimplexFlagNull );

        M.check_dirichlet_flags();

        
        LOG << "Prepare scalar fields for testing..." << nl;
        

        auto trafo     = [](const FloatVector& vec) -> FloatVector { 
            assert( vec.getdimension() == 3 );
            // return vec; // TODO: correct trafo
            return vec.sumnorm() / vec.l2norm() * vec;
        };

        auto trafo_inv = [](const FloatVector& vec) -> FloatVector { 
            assert( vec.getdimension() == 3 );
            // return vec; // TODO: correct trafo
            return vec.l2norm() / vec.sumnorm() * vec;
        };

        // auto trafo_jacobian
        auto jacobian = [](const FloatVector& vec) -> DenseMatrix { 
            assert( vec.getdimension() == 3 );
            
            // return scalar * DenseMatrix(3,3, kronecker<int> ); // TODO: deactivate 

            Float x = vec[0];
            Float y = vec[1];
            Float z = vec[2];

            Float ax = absolute(x);
            Float ay = absolute(y);
            Float az = absolute(z);

            Float sx = sign(x);
            Float sy = sign(y);
            Float sz = sign(z);

            Float dx = ( ax >= ay and ax >= az ) ? sx : 0.;
            Float dy = ( ay >= ax and ay >= az ) ? sy : 0.;
            Float dz = ( az >= ax and az >= ay ) ? sz : 0.;
            
            Float li  = vec.sumnorm();
            Float l2  = vec.l2norm();
            Float l2c = l2*l2*l2;

            Float liol2  = li / l2;
            Float liol2c = li / l2c;
            
            return DenseMatrix( 3, 3, {
                liol2 + x * dx/l2 - x*x* liol2c,    x * dy/l2 - x * y * liol2c,        x * dz/l2 - x * z * liol2c, 
                y * dx/l2 - y * x * liol2c,         liol2 + x * dx/l2 - x*x* liol2c,   y * dz/l2 - y * z * liol2c,  
                z * dx/l2 - z * z * liol2c,         z * dy/l2 - z * y * liol2c,        liol2 + z * dz/l2 - z*z* liol2c, 
            });
        };

        // jacobian matrix of { x * max(x,y,z) / sqrt(xx+yy+zz), y * max(x,y,z) / sqrt(xx+yy+zz), z * max(x,y,z) / sqrt(xx+yy+zz) }


        
        auto physical_f = 
            [=](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ 
                    1.
                    });
            };
        
        auto parametric_f = 
            [=](const FloatVector& vec) -> FloatVector{
                return physical_f( trafo(vec) );
            };
        
        
        auto physical_A = 
            [=](const FloatVector& vec) -> DenseMatrix{
                assert( vec.getdimension() == 3 );
                Float scalar = vec.l2norm() > 0.5 ? 1. : 2. ;
                return scalar * DenseMatrix(3,3, kronecker<int> );
            };
        
        auto parametric_A = 
            [=](const FloatVector& vec) -> DenseMatrix{
                return physical_A( trafo(vec) );
            };
        
        // std::function<DenseMatrix(const FloatVector&)> 
        auto weight_scalar = 
            [=](const FloatVector& vec) -> DenseMatrix{
                assert( vec.getdimension() == 3 );
                // return DenseMatrix(1,1, kronecker<int> );
                return DenseMatrix(1,1,absolute(Determinant(jacobian(vec))));
            };
        
        
        // std::function<DenseMatrix(const FloatVector&)> 
        auto weight_vector = 
            [=](const FloatVector& vec) -> DenseMatrix{
                assert( vec.getdimension() == 3 );
                // return DenseMatrix(3,3, kronecker<int> );
                auto jac = jacobian(vec);

                auto det = Determinant(jac);

                auto invjac = Inverse( jac );

                return absolute(det) * invjac * parametric_A(vec) * Transpose(invjac);
            };
        
        
        
        

        

        LOG << "Solving Poisson Problem with Neumann boundary conditions" << nl;

        const int min_l = 1; 
        const int max_l = 5;

        const int min_r = 1;
        const int max_r = 1;

        const int r_plus = 2;
        const int w_plus = 2;
        
        ConvergenceTable contable("Mass error");
        
        contable << "u_error" << "du_error" << "residual" << "time" << nl;
        

        assert( 0 <= min_l and min_l <= max_l );
        assert( 0 <= min_r and min_r <= max_r );
        
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ ){
            
            LOG << "Level: " << l << "/" << max_l << nl;
            LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            // if( l != 0 )
            for( int r = min_r; r <= max_r; r++ ) 
            {
                
                int w = r;

                LOG << "...assemble scalar mass matrices" << nl;
        
                SparseMatrix     scalar_massmatrix = FEECBrokenCoefficientMassMatrix( M, M.getinnerdimension(), 0, r         , w         , weight_scalar );
                SparseMatrix aug_scalar_massmatrix = FEECBrokenCoefficientMassMatrix( M, M.getinnerdimension(), 0, r + r_plus, w + w_plus, weight_scalar );

                LOG << "...assemble vector mass matrix" << nl;
        
                SparseMatrix     vector_massmatrix = FEECBrokenCoefficientMassMatrix( M, M.getinnerdimension(), 1, r-1         , w         , weight_vector );
                SparseMatrix aug_vector_massmatrix = FEECBrokenCoefficientMassMatrix( M, M.getinnerdimension(), 1, r-1 + r_plus, w + w_plus, weight_vector );
                
                LOG << "...assemble differential matrix and transpose" << nl;

                SparseMatrix     diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r          );
                SparseMatrix aug_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r + r_plus );

                SparseMatrix     diffmatrix_t =     diffmatrix.getTranspose();
                SparseMatrix aug_diffmatrix_t = aug_diffmatrix.getTranspose();

                LOG << "...assemble inclusion matrix and transpose" << nl;
        
                SparseMatrix     incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r          );
                SparseMatrix aug_incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r + r_plus );
                
                SparseMatrix     incmatrix_t =     incmatrix.getTranspose();
                SparseMatrix aug_incmatrix_t = aug_incmatrix.getTranspose();

                LOG << "...assemble elevation matrix" << nl;
        
                SparseMatrix elevation_matrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_plus );

                LOG << "...assemble stiffness matrix" << nl;
        
                auto opr  = diffmatrix & incmatrix;
                auto opl  = opr.getTranspose(); 
                auto stiffness = opl & ( vector_massmatrix & opr );                
                stiffness.sortentries();
                auto stiffness_csr = MatrixCSR( stiffness );
                
                auto aug_opr  = aug_diffmatrix & aug_incmatrix;
                auto aug_opl  = aug_opr.getTranspose(); 
                auto aug_stiffness = aug_opl & ( aug_vector_massmatrix & aug_opr );                
                aug_stiffness.sortentries();
                auto aug_stiffness_csr = MatrixCSR( aug_stiffness );
                
                {

                    const auto& function_rhs  = parametric_f;
                    
                    LOG << "...interpolate explicit solution and rhs" << nl;
        
                    FloatVector     interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,        function_rhs  );
                    FloatVector aug_interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r+r_plus, function_rhs  );

                    FloatVector     rhs =     incmatrix_t * (     scalar_massmatrix *     interpol_rhs );
                    FloatVector aug_rhs = aug_incmatrix_t * ( aug_scalar_massmatrix * aug_interpol_rhs );

                    FloatVector     sol(     incmatrix.getdimin(), 0. );
                    FloatVector aug_sol( aug_incmatrix.getdimin(), 0. );
                    
                

                    LOG << "...iterative solver" << nl;
                    
                    timestamp start = gettimestamp();
                    
                    assert( stiffness_csr.getC() );


                    {
                        sol.zero();

                        auto diagonal = stiffness.diagonal();
                        FloatVector residual( rhs );
                        ConjugateGradientSolverCSR_SSOR( 
                            sol.getdimension(), 
                            sol.raw(), 
                            rhs.raw(), 
                            stiffness_csr.getA(), stiffness_csr.getC(), stiffness_csr.getV(),
                            residual.raw(),
                            desired_precision,
                            0,
                            diagonal.raw(),
                            1.0
                        );

                    }

                    {
                        aug_sol.zero();

                        auto aug_diagonal = aug_stiffness.diagonal();
                        FloatVector aug_residual( aug_rhs );
                        ConjugateGradientSolverCSR_SSOR( 
                            aug_sol.getdimension(), 
                            aug_sol.raw(), 
                            aug_rhs.raw(), 
                            aug_stiffness_csr.getA(), aug_stiffness_csr.getC(), aug_stiffness_csr.getV(),
                            aug_residual.raw(),
                            desired_precision,
                            0,
                            aug_diagonal.raw(),
                            1.0
                        );

                    }

                    timestamp end = gettimestamp();
                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                
                    
                    
                    LOG << "...compute error and residual:" << nl;

                    FloatVector error     = aug_incmatrix * aug_sol - elevation_matrix * incmatrix * sol;
                    FloatVector graderror = aug_diffmatrix * ( aug_incmatrix * aug_sol - elevation_matrix * incmatrix * sol );
                    Float errornorm       = std::sqrt( error * ( aug_scalar_massmatrix * error ) );
                    Float graderrornorm   = std::sqrt( graderror * ( aug_vector_massmatrix * graderror ) );
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
                        
                        auto computed_grad = diffmatrix * incmatrix * sol;
                        
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
