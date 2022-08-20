

/**/

#include <fstream>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
    
    LOG << "Unit Test for Solution of Darcy Problem" << endl;
    
    // LOG << std::setprecision(10);

    if(true){

        LOG << "Initial mesh..." << endl;
        
        MeshSimplicial2D M = StandardSquare2D_simple();
        
        M.check();
        
        LOG << "Prepare scalar fields for testing..." << endl;
        

        std::function<FloatVector(const FloatVector&)> constant_one
            = [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 1. });
                };
        
        
        
        


        
        // std::function<FloatVector(const std::function<FloatVector(const FloatVector&) ) >scalarfield = 
        
        const Float xfeq = 1.;
        const Float yfeq = 1.;
        
        
        // u dx + v dy -> u_y dydx + v_x dxdy = ( v_x - u_y ) dxdy
        
        // phi -> ( - phi_y, phi_x ) -> ( - phi_xx - phi_yy ) dxdy
        
        std::function<FloatVector(const FloatVector&)> experiment_sol = 
            [=](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                // return FloatVector({ 1. });
                return FloatVector({ std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) });
            };
        
        
        std::function<FloatVector(const FloatVector&)> experiment_grad = 
            [=](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                // return FloatVector({ 1. });
                return FloatVector( { 
                            yfeq * Constants::twopi * std::sin( xfeq * Constants::twopi * vec[0] ) * std::cos( yfeq * Constants::twopi * vec[1] ), 
                        -xfeq * Constants::twopi * std::cos( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ),
                    });
            };
        
        
        std::function<FloatVector(const FloatVector&)> experiment_rhs = 
            [=](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ 
                    xfeq*xfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] )
                    +
                    yfeq*yfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] )
                    });
            };
        
        
        

        

        LOG << "Solving Poisson Problem with Neumann boundary conditions" << endl;

        const int min_l = 0; 
        const int max_l = 3;
        
        const int min_r = 5;
        const int max_r = 5;

        const int aug_r = 1;
        
        
        ConvergenceTable contable("Mass error");
        
        contable << "sigma_error" << "u_error" << "residual" << "time";
        

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
                
                LOG << "... assemble matrices" << endl;
        
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
                
                SparseMatrix vector_massmatrix_inv = FEECBrokenMassMatrix_cellwiseinverse( M, M.getinnerdimension(), 1, r   );
                
                SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

                SparseMatrix diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
                SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                SparseMatrix volume_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r-1 );
                SparseMatrix volume_incmatrix_t = volume_incmatrix.getTranspose();

                auto mat_A  = vector_incmatrix_t & vector_massmatrix & vector_incmatrix;
                mat_A.sortandcompressentries();
                
                auto mat_Bt = vector_incmatrix_t & diffmatrix_t & volume_massmatrix & volume_incmatrix; // upper right
                mat_Bt.sortandcompressentries();
                
                auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower bottom
                mat_B.sortandcompressentries();
                
                auto A  = MatrixCSR( mat_A  );
                auto Bt = MatrixCSR( mat_Bt );
                auto B  = MatrixCSR( mat_B  );
                
                auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
                
                auto Schur = B * inv(A,1e-14) * Bt;

                auto PA = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )
                              + MatrixCSR( vector_incmatrix_t & diffmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix );
                auto PC = MatrixCSR( volume_incmatrix_t & volume_massmatrix & volume_incmatrix );
                    
                
                {

                    const auto& function_sol  = experiment_sol;
                    const auto& function_grad = experiment_grad;
                    const auto& function_rhs  = experiment_rhs;
                    
                    LOG << "...interpolate explicit solution and rhs" << endl;
                    
                    FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r,   function_grad );
                    FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 2, r-1, function_sol  );
                    FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 2, r-1, function_rhs  );
                    
                    
                    {
                        
                        FloatVector rhs = volume_incmatrix_t * ( volume_massmatrix * interpol_rhs );

                        FloatVector sol( volume_incmatrix.getdimin(), 0. );
                        
                        LOG << "...iterative solver" << endl;
                        
                        
                        sol.zero();
                        
                        FloatVector res = sol;
                        


                        timestamp start = gettimestamp();

                        {

                            const auto PAinv = inv(PA,desired_precision,-1);
                            const auto PCinv = inv(PC,desired_precision,-1);

                            FloatVector  x_A( A.getdimin(),  0. ); 
                            FloatVector& x_C = sol;
                            
                            const FloatVector  b_A( A.getdimin(),  0. ); 
                            const FloatVector& b_C = rhs; 
                            
                            auto Z  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix

                            BlockHerzogSoodhalterMethod( 
                                x_A, 
                                x_C, 
                                b_A, 
                                b_C, 
                                -A, Bt, B, Z, 
                                desired_precision,
                                1,
                                PAinv, PCinv
                            );

                        }

                        timestamp end = gettimestamp();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                        
                        
                        auto grad = inv(A,desired_precision) * Bt * sol;

                        LOG << "...compute error and residual:" << endl;

                        
                        FloatVector interpol_grad_aug = Interpolation( M, M.getinnerdimension(), 1, r + aug_r,     function_grad );
                        FloatVector interpol_sol_aug  = Interpolation( M, M.getinnerdimension(), 2, r + aug_r - 1, function_sol  );

                        SparseMatrix vector_elevation_matrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r,   aug_r );
                        SparseMatrix volume_elevation_matrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r-1, aug_r );

                        SparseMatrix vector_massmatrix_aug = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r + aug_r     );
                        SparseMatrix volume_massmatrix_aug = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r + aug_r - 1 );

                        auto errornorm_aux_sol  = interpol_sol_aug  - volume_elevation_matrix * volume_incmatrix *  sol;
                        auto errornorm_aux_grad = interpol_grad_aug - vector_elevation_matrix * vector_incmatrix * grad;

                        Float errornorm_sol  = sqrt( errornorm_aux_sol  * ( volume_massmatrix_aug *  errornorm_aux_sol ) );
                        Float errornorm_grad = sqrt( errornorm_aux_grad * ( vector_massmatrix_aug * errornorm_aux_grad ) );
                        Float residualnorm   = ( rhs - B * inv(A,1e-14) * Bt * sol ).norm();

                        LOG << "error:     " << errornorm_sol << endl;
                        LOG << "aux error: " << errornorm_grad << endl;
                        LOG << "residual:  " << residualnorm << endl;

                        contable << errornorm_sol;
                        contable << errornorm_grad;
                        contable << residualnorm;
                        contable << Float( end - start );
                        contable << nl;

                    }





                    contable.lg();
    
                }
                
            }

            if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
            
            

        } 
    
    }
    
    
    
    
    LOG << "Finished Unit Test" << endl;
    
    return 0;
}
