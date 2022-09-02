

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
#include "../../fem/local.polynomialmassmatrix.hpp"
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
        
        
        ConvergenceTable contable("Mass error");
        
        contable << "sigma_error" << "u_error" << "residual" << "time";
        

        assert( 0 <= min_l and min_l <= max_l );
        assert( 0 <= min_r and min_r <= max_r );
            
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ ){
            
            LOG << "Level: " << l << "/" << max_l << nl;
            LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            if( l != 0 )
            for( int r = min_r; r <= max_r; r++ ) 
            {
                
                LOG << "... assemble matrices" << nl; // TODO: correct the degrees, perhaps via degree elevation
        
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
                
//                 SparseMatrix vector_massmatrix_inv = FEECBrokenMassMatrix_cellwiseinverse( M, M.getinnerdimension(), 1, r );
                
                SparseMatrix volume_elevationmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r-1, 1 );
                SparseMatrix volume_elevationmatrix_t = volume_elevationmatrix.getTranspose();

                SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );

                SparseMatrix diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                SparseMatrix vector_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r );
                SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                SparseMatrix volume_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 2, r );
                SparseMatrix volume_incmatrix_t = volume_incmatrix.getTranspose();

                auto mat_A  = vector_incmatrix_t & vector_massmatrix & vector_incmatrix;
                mat_A.sortandcompressentries();
                
                auto mat_Bt = vector_incmatrix_t & diffmatrix_t & volume_elevationmatrix_t & volume_massmatrix & volume_incmatrix; // upper right
                mat_Bt.sortandcompressentries();
                
                auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower bottom
                mat_B.sortandcompressentries();
                
                auto A  = MatrixCSR( mat_A  );
                auto Bt = MatrixCSR( mat_Bt );
                auto B  = MatrixCSR( mat_B  );
                
                auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
                
                auto Schur = B * inv(A,1e-10) * Bt;
                
                {

                    const auto& function_sol  = experiment_sol;
                    const auto& function_grad = experiment_grad;
                    const auto& function_rhs  = experiment_rhs;
                    
                    LOG << "...interpolate explicit solution and rhs" << nl;
                    
                    FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r, function_grad );
                    FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 2, r-1, function_sol  );
                    FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 2, r-1, function_rhs  );
                    
                    timestamp start = gettimestamp();

                    FloatVector sol( volume_incmatrix.getdimin(), 0. );
                    sol.zero();

                    FloatVector rhs = volume_incmatrix_t * ( volume_massmatrix * volume_elevationmatrix * interpol_rhs );

                    {
                        LOG << "...iterative solver" << nl;

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
                        
                    }

                    timestamp end = gettimestamp();
                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                    LOG << "...compute error and residual:" << nl;

                    auto grad = inv(A,1e-14) * Bt * sol;

                    auto errornorm_aux_sol  = volume_elevationmatrix * interpol_sol  - volume_incmatrix * sol;
                    auto errornorm_aux_grad = interpol_grad - vector_incmatrix * grad;

                    Float errornorm_sol  = sqrt( errornorm_aux_sol  * ( volume_massmatrix *  errornorm_aux_sol ) );
                    Float errornorm_grad = sqrt( errornorm_aux_grad * ( vector_massmatrix * errornorm_aux_grad ) );
                    Float residualnorm   = ( rhs - B * inv(A,1e-10) * Bt * sol ).norm();

                    LOG << "error:     " << errornorm_sol << nl;
                    LOG << "aux error: " << errornorm_grad << nl;
                    LOG << "residual:  " << residualnorm << nl;

                    contable << errornorm_sol;
                    contable << errornorm_grad;
                    contable << residualnorm;
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
