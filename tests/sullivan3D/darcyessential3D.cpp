

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../solver/iterativesolver.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/minres.hpp"
// #include "../../solver/herzogsoodhalter.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
    
    LOG << "Unit Test for Solution of Darcy Problem" << nl;
    
    // LOG << std::setprecision(10);

    if(true){

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial3D M = StandardCube3D();
        
        M.check();
        
        M.automatic_dirichlet_flags();
        M.check_dirichlet_flags();

        
        LOG << "Prepare scalar fields for testing..." << nl;
        

        // std::function<FloatVector(const std::function<FloatVector(const FloatVector&) ) >scalarfield = 
        
        const Float xfeq = 1.;
        const Float yfeq = 2.;
        const Float zfeq = 3.;
        
        
        // u dx + v dy -> u_y dydx + v_x dxdy = ( v_x - u_y ) dxdy
        
        // phi -> ( - phi_y, phi_x ) -> ( - phi_xx - phi_yy ) dxdy
        
        std::function<FloatVector(const FloatVector&)> experiment_sol = 
            [=](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] ) });
                // return FloatVector({ blob( vec[0] ) * blob( vec[1] ) * blob( vec[2] ) });
            };
        
        
        std::function<FloatVector(const FloatVector&)> experiment_grad = 
            [=](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector( { 
                        -zfeq * Constants::twopi * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::cos( zfeq * Constants::twopi * vec[2] ), //xy
                        +yfeq * Constants::twopi * std::sin( xfeq * Constants::twopi * vec[0] ) * std::cos( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] ), //xz
                        -xfeq * Constants::twopi * std::cos( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] )  //yz
                        // - blob( vec[0] )     * blob( vec[1] )     * blob_dev( vec[2] ),
                        // + blob( vec[0] )     * blob_dev( vec[1] ) * blob( vec[2] ),
                        // - blob_dev( vec[0] ) * blob( vec[1] )     * blob( vec[2] )
                    });
            };
        
        
        std::function<FloatVector(const FloatVector&)> experiment_rhs = 
            [=](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ 
                    xfeq*xfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] )
                    +
                    yfeq*yfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] )
                    +
                    zfeq*zfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] )
                    // -
                    //  blob_devdev( vec[0] ) * blob( vec[1] ) * blob( vec[2] )
                    // -
                    // blob( vec[0] ) * blob_devdev( vec[1] ) * blob( vec[2] )
                    // -
                    // blob( vec[0] ) * blob( vec[1] ) * blob_devdev( vec[2] )
                    });
            };
        
        
        

        

        const int min_l = 0; 
        const int max_l = 5;
        
        const int min_r = 1;
        const int max_r = 1;

        const int r_plus_pseudo = 0;
        const int r_plus_volume = 0;
        
        
        ConvergenceTable contable("Mass error");
        
        contable << "sigma_error" << "u_error" << "sigma_res" << "u_res" << "time";
        

        assert( 0 <= min_l and min_l <= max_l );
        assert( 0 <= min_r and min_r <= max_r );
        assert( 0 <= r_plus_pseudo and 0 <= r_plus_volume );
            
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ ){
            
            LOG << "Level: " << l << "/" << max_l << nl;
            LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            if( l != 0 )
            for( int r = min_r; r <= max_r; r++ ) 
            {
                
                LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                    
                LOG << "integration with: " << r_plus_pseudo << ", " << r_plus_volume << nl;
                    
                LOG << "... assemble matrices" << nl;
        
                SparseMatrix pseudo_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r   + r_plus_pseudo );
                
                SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r-1 + r_plus_volume );

                SparseMatrix diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 2, r );
                SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                SparseMatrix pseudo_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r   );
                SparseMatrix pseudo_incmatrix_t = pseudo_incmatrix.getTranspose();

                SparseMatrix volume_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 3, r-1 );
                SparseMatrix volume_incmatrix_t = volume_incmatrix.getTranspose();

                SparseMatrix pseudo_elevmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r  , r_plus_pseudo );
                SparseMatrix pseudo_elevmatrix_t = pseudo_elevmatrix.getTranspose();

                SparseMatrix volume_elevmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 3, r-1, r_plus_volume );
                SparseMatrix volume_elevmatrix_t = volume_elevmatrix.getTranspose();

                auto mat_A  = pseudo_incmatrix_t & pseudo_elevmatrix_t & pseudo_massmatrix & pseudo_elevmatrix & pseudo_incmatrix;
                mat_A.sortandcompressentries();
                
                auto mat_Bt = pseudo_incmatrix_t & diffmatrix_t & volume_elevmatrix_t & volume_massmatrix & volume_elevmatrix & volume_incmatrix; // upper right
                mat_Bt.sortandcompressentries();
                
                auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & pseudo_incmatrix; // lower bottom
                mat_B.sortandcompressentries();
                
                auto A  = MatrixCSR( mat_A  );
                auto Bt = MatrixCSR( mat_Bt );
                auto B  = MatrixCSR( mat_B  );
                
                auto negA = - A; 
                
                // auto Schur = B * inv(A,1e-10) * Bt;
                
                {

                    const auto& function_sol  = experiment_sol;
                    const auto& function_grad = experiment_grad;
                    const auto& function_rhs  = experiment_rhs;
                    
                    LOG << "...interpolate explicit solution and rhs" << nl;
                    
                    FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 2, r   + r_plus_pseudo, function_grad );
                    FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 3, r-1 + r_plus_volume, function_sol  );
                    FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 3, r-1 + r_plus_volume, function_rhs  );
                    
                    FloatVector rhs = volume_incmatrix_t * volume_elevmatrix_t * ( volume_massmatrix * interpol_rhs );

                    FloatVector sol( volume_incmatrix.getdimin(), 0. );
                    
                    timestamp start = gettimestamp();

                    //{

                        LOG << "...iterative solver" << nl;
                        
                        auto PA = MatrixCSR( pseudo_incmatrix_t & pseudo_elevmatrix_t & pseudo_massmatrix & pseudo_elevmatrix & pseudo_incmatrix )
                                + MatrixCSR( pseudo_incmatrix_t & diffmatrix_t & volume_elevmatrix_t & volume_massmatrix & volume_elevmatrix & diffmatrix & pseudo_incmatrix );
                        auto PC = MatrixCSR( volume_incmatrix_t & volume_elevmatrix_t & volume_massmatrix & volume_elevmatrix & volume_incmatrix );
                
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
                            negA, Bt, B, Z, 
                            desired_precision,
                            1,
                            PAinv, PCinv
                        );

                    //}
                    
                    timestamp end = gettimestamp();
                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                                        
                    auto grad = x_A; // inv(A,desired_precision) * Bt * sol;

                    LOG << "...compute error and residual:" << nl;

                        
                    // improved error estimation 
                    
                    FloatVector interpol_grad_aug = Interpolation( M, M.getinnerdimension(), 2, r   + r_plus_pseudo, function_grad );
                    FloatVector interpol_sol_aug  = Interpolation( M, M.getinnerdimension(), 3, r-1 + r_plus_volume, function_sol  );

                    auto errornorm_aux_grad = interpol_grad_aug - pseudo_elevmatrix * pseudo_incmatrix * grad;
                    auto errornorm_aux_sol  = interpol_sol_aug  - volume_elevmatrix * volume_incmatrix *  sol;

                    Float errornorm_grad = sqrt( errornorm_aux_grad * ( pseudo_massmatrix * errornorm_aux_grad ) );
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
                    
                }
                
                contable.lg();
                
            }

            if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
    
        } 
    
    }
    
    LOG << "Finished Unit Test" << nl;
    
    return 0;
}
