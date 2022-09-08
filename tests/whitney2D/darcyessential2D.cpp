

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../solver/iterativesolver.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/minres.hpp"
// #include "../../solver/herzogsoodhalter.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
    
    LOG << "Unit Test for Solution of Darcy Problem" << nl;
    
    // LOG << std::setprecision(10);

    if(true){

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
//             M.automatic_dirichlet_flags();
//             M.check_dirichlet_flags();

        
        LOG << "Prepare scalar fields for testing..." << nl;
        

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
        
        
        

        

        LOG << "Solving Poisson Problem with Neumann boundary conditions" << nl;

        const int min_l = 0; 
        const int max_l = 5;
        
        const int min_r = 1;
        const int max_r = 1;
        
        
        const int r_plus_vector = 1;
        const int r_plus_volume = 1;
        
        
        ConvergenceTable contable("Mass error");
        
        contable << "sigma_error" << "u_error" << "sigma_res" << "u_res" << "time";
        

        assert( 0 <= min_l and min_l <= max_l );
        assert( 0 <= min_r and min_r <= max_r );
        assert( 0 <= r_plus_vector and 0 <= r_plus_volume );
            
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ ){
            
            LOG << "Level: " << l << "/" << max_l << nl;
            LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            if( l != 0 )
            for( int r = min_r; r <= max_r; r++ ) 
            {
                
                LOG << "... assemble matrices" << nl; // TODO: correct the degrees, perhaps via degree elevation
        
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r + r_plus_vector );
                
                SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r + r_plus_volume );

                SparseMatrix vector_augdiffmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r-1, 1 )
                                                      &
                                                      FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r )
                                                      ;
                SparseMatrix vector_augdiffmatrix_t = vector_augdiffmatrix.getTranspose();

                SparseMatrix vector_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r );
                SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                SparseMatrix volume_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 2, r );
                SparseMatrix volume_incmatrix_t = volume_incmatrix.getTranspose();

                SparseMatrix vector_elevmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r, r_plus_vector );
                SparseMatrix vector_elevmatrix_t = vector_elevmatrix.getTranspose();

                SparseMatrix volume_elevmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r, r_plus_volume );
                SparseMatrix volume_elevmatrix_t = volume_elevmatrix.getTranspose();

                auto mat_A  = vector_incmatrix_t & vector_elevmatrix_t & vector_massmatrix & vector_elevmatrix & vector_incmatrix;
                mat_A.sortandcompressentries();
                
                auto mat_Bt = vector_incmatrix_t & vector_augdiffmatrix_t & volume_elevmatrix_t & volume_massmatrix & volume_elevmatrix & volume_incmatrix; // upper right
                mat_Bt.sortandcompressentries();
                
                auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix; // lower bottom
                mat_B.sortandcompressentries();
                
                auto A  = MatrixCSR( mat_A  );
                auto Bt = MatrixCSR( mat_Bt );
                auto B  = MatrixCSR( mat_B  );
                
                auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
                
                {

                    const auto& function_sol  = experiment_sol;
                    const auto& function_grad = experiment_grad;
                    const auto& function_rhs  = experiment_rhs;
                    
                    LOG << "...interpolate explicit solution and rhs" << nl;
                    
                    FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r + r_plus_vector, function_grad );
                    FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 2, r + r_plus_volume, function_sol  );
                    FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 2, r + r_plus_volume, function_rhs  );
                    
                    timestamp start = gettimestamp();

                    FloatVector sol( volume_incmatrix.getdimin(), 0. );
                    sol.zero();

                    FloatVector rhs = volume_incmatrix_t * volume_elevmatrix_t 
                        * ( volume_massmatrix * interpol_rhs );

                    // {
                        LOG << "...iterative solver" << nl;

                        // auto PA = MatrixCSR( vector_incmatrix_t & vector_elevmatrix_t & vector_massmatrix & vector_elevmatrix & vector_incmatrix )
                        //         + MatrixCSR( vector_incmatrix_t & vector_augdiffmatrix_t & volume_elevmatrix_t & volume_massmatrix & volume_elevmatrix & vector_augdiffmatrix & vector_incmatrix );
                        // auto PC = MatrixCSR( volume_incmatrix_t & volume_elevmatrix_t & volume_massmatrix & volume_elevmatrix & volume_incmatrix );


                        FloatVector res = sol;
                        
                        HodgeConjugateResidualSolverCSR_SSOR( 
                            B.getdimout(), 
                            A.getdimout(), 
                            sol.raw(), 
                            rhs.raw(), 
                            A.getA(),   A.getC(),  A.getV(), 
                            B.getA(),   B.getC(),  B.getV(), 
                            Bt.getA(), Bt.getC(), Bt.getV(), 
                            C.getA(),   C.getC(),  C.getV(), 
                            res.raw(),
                            1e-10,
                            1,
                            desired_precision,
                            0
                        );
                        
                    // }

                    timestamp end = gettimestamp();
                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                    auto grad = inv(A,desired_precision) * Bt * sol;

                    LOG << "...compute error and residual:" << nl;

                        
                    // improved error estimation 
                    
                    FloatVector interpol_grad_aug = Interpolation( M, M.getinnerdimension(), 1, r + r_plus_vector, function_grad );
                    FloatVector interpol_sol_aug  = Interpolation( M, M.getinnerdimension(), 2, r + r_plus_volume, function_sol  );

                    auto errornorm_aux_grad = interpol_grad_aug - vector_elevmatrix * vector_incmatrix * grad;
                    auto errornorm_aux_sol  = interpol_sol_aug  - volume_elevmatrix * volume_incmatrix *  sol;

                    Float errornorm_grad = sqrt( errornorm_aux_grad * ( vector_massmatrix * errornorm_aux_grad ) );
                    Float errornorm_sol  = sqrt( errornorm_aux_sol  * ( volume_massmatrix *  errornorm_aux_sol ) );
                    
                    Float residual_sol   = ( rhs - B * grad ).norm();
                    Float residual_grad  = ( - A * grad + Bt * sol ).norm();

                    LOG << "error:     " << errornorm_sol << nl;
                    LOG << "aux error: " << errornorm_grad << nl;
                    LOG << "residual:  " << residual_sol << nl;
                    LOG << "residual:  " << residual_grad << nl;

                    contable << errornorm_sol;
                    contable << errornorm_grad;
                    contable << residual_sol;
                    contable << residual_grad;
                    contable << Float( end - start );
                    contable << nl;
                    
                    contable.lg();
                    
                }
                
            }

            if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
            
        } 
    
    }
    
    LOG << "Finished Unit Test" << nl;
    
    return 0;
}
