

/**/

#include <ostream>
#include <fstream>
// #include <iomanip>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
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
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        M.automatic_dirichlet_flags();
        M.check_dirichlet_flags();

        
        LOG << "Prepare scalar fields for testing..." << endl;
        

        std::function<FloatVector(const FloatVector&)> constant_one
            = [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 1. });
                };
        
        
        // u dx + v dy -> u_y dydx + v_x dxdy = - u_y dxdy + v_x dxdy = ( v_x - u_y ) dxdy
        // 
        // \int ( v_x - u_y ) phi dxdy
        //    = \int ( v_x - u_y ) phi 
        //    = \int -v phi_x + u phi_y
        //    = \int u phi_y - v phi_x
        //    = \int (u,v) * ( phi_y, -phi_x )
        //    = \int u dx phi_y dy - v dy phi_x dx
        //    = \int u phi_y dxdy - v phi_x dy dx
        //    = \int u phi_y dxdy + v phi_x dxdy
        //      
        //      phi -> ( phi_y, -phi_x )
        // 
        // TOT = ( v_xy - u_yy, - v_xx + u_xy ) 
        // TOT = ( - u_yy + v_xy , - v_xx + u_xy ) 
        //      
        // phi -> ( - phi_y, phi_x ) -> ( - phi_xx - phi_yy ) dxdy
        
        // u dx + v dy -> ( v_x - u_y ) dxdy -> ( - v_xy + u_yy, v_xx - u_xy )
        
        std::function<FloatVector(const FloatVector&)> experiment_sol = 
            [=](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                // return FloatVector({ 1. });
                return FloatVector({ 
                    0 //bumpfunction_dev(vec[0])*bumpfunction(vec[1])
                    ,
                    0 // bumpfunction(vec[0])*bumpfunction_dev(vec[1])
                    });
            };
        
        std::function<FloatVector(const FloatVector&)> experiment_rhs = 
            [=](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ 
                    bumpfunction_dev(vec[0])*bumpfunction(vec[1])
                    ,
                    bumpfunction(vec[0])*bumpfunction_dev(vec[1])
                    // - bumpfunction(vec[0])*bumpfunction_devdev(vec[1]) + bumpfunction_dev(vec[0])*bumpfunction_dev(vec[1])
                    // ,
                    // - bumpfunction_devdev(vec[0])*bumpfunction(vec[1]) + bumpfunction_dev(vec[0])*bumpfunction_dev(vec[1])
                    });
            };
        
        std::function<FloatVector(const FloatVector&)> experiment_aux = 
            [=](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                // return FloatVector({ 1. });
                return FloatVector( { 
                    bumpfunction_dev(vec[0])*bumpfunction(vec[1])
                    });
            };
        
        
        
        
        

        

        LOG << "Solving Poisson Problem with Neumann boundary conditions" << endl;

        const int min_l = 0; 
        const int max_l = 5;
        
        const int min_r = 2;
        const int max_r = 2;
        
        
        ConvergenceTable contable("Mass error and numerical residuals");
        
        contable << "u_error" << "sigma_error" << "u_res" << "sigma_res";
        

        assert( 0 <= min_l and min_l <= max_l );
        assert( 0 <= min_r and min_r <= max_r );
            
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ ){
            
            LOG << "Level: " << l << "/" << max_l << std::endl;
            LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            for( int r = min_r; r <= max_r; r++ ) 
            {
                
                LOG << "... assemble matrices" << endl;
        
                SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 );
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
                SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

                SparseMatrix scalar_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 );
                SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

                SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 );
                SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
                SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

                auto mat_A  = vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix;
                mat_A.sortandcompressentries();
                    
                auto mat_Bt = vector_incmatrix_t & vector_massmatrix & scalar_diffmatrix & scalar_incmatrix; // upper right
                mat_Bt.sortandcompressentries();
                
                auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower bottom
                mat_B.sortandcompressentries();
                
                auto A  = MatrixCSR( mat_A  );
                auto Bt = MatrixCSR( mat_Bt );
                auto B  = MatrixCSR( mat_B  );
                
                auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
                
                
                const auto& function_sol = experiment_sol;
                const auto& function_rhs = experiment_rhs;
                const auto& function_aux = experiment_aux;
                
                LOG << "...interpolate explicit solution and rhs" << endl;
                
                FloatVector interpol_sol = Interpolation( M, M.getinnerdimension(), 1, r,   function_sol );
                FloatVector interpol_rhs = Interpolation( M, M.getinnerdimension(), 1, r,   function_rhs );
                FloatVector interpol_aux = Interpolation( M, M.getinnerdimension(), 0, r+1, function_aux );
                
                FloatVector rhs_sol = vector_incmatrix_t * vector_massmatrix * interpol_rhs;
                FloatVector rhs_aux = scalar_incmatrix_t * scalar_diffmatrix_t * vector_massmatrix * interpol_sol ;// FloatVector( B.getdimout(), 0. );

                LOG << "...measure interpolation commutativity" << endl;
    
                auto X = Block2x2Operator( A.getdimout() + B.getdimout(), A.getdimin() + Bt.getdimin(), A, Bt, B, C );

                FloatVector sol_whole( A.getdimin()  + Bt.getdimin(),  0. );
                FloatVector rhs_whole( A.getdimout() +  B.getdimout(), 0. );
                
                sol_whole.zero();
                rhs_whole.setslice(             0, rhs_sol );
                rhs_whole.setslice( A.getdimout(), rhs_aux );
                
                HerzogSoodhalterMethod Solver( X );
                Solver.print_modulo        = 500;
                Solver.max_iteration_count = 2 * sol_whole.getdimension();

                timestamp start = gettimestamp();
                Solver.solve( sol_whole, rhs_whole );
                timestamp end = gettimestamp();

                LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                LOG << "...compute error and residual:" << endl;

                FloatVector res_whole = ( rhs_whole - X * sol_whole );

                Float residualnorm = res_whole.norm();

                LOG << "combined system residual:  " << residualnorm << endl;

                FloatVector sol = sol_whole.getslice(             0, A.getdimout() );
                FloatVector aux = sol_whole.getslice( A.getdimout(), B.getdimout() );
                
                LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                
                
                LOG << "...compute error and residual:" << endl;

                auto errornorm_aux_sol = interpol_sol - vector_incmatrix * sol;
                auto errornorm_aux_aux = interpol_aux - scalar_incmatrix * aux;

                Float errornorm_sol = sqrt( errornorm_aux_sol * ( vector_massmatrix * errornorm_aux_sol ) );
                Float errornorm_aux = sqrt( errornorm_aux_aux * ( scalar_massmatrix * errornorm_aux_aux ) );
                
                Float residual_sol  = ( rhs_sol - A * sol - Bt * aux ).norm();
                Float residual_aux  = ( rhs_aux - B * sol -  C * aux ).norm();

                LOG << "error:     " << errornorm_sol << endl;
                LOG << "aux error: " << errornorm_aux << endl;
                LOG << "residual:  " << residualnorm << endl;

                contable << errornorm_sol;
                contable << errornorm_aux;
                
                contable << residual_sol;
                contable << residual_aux;
                
                contable << nl;

                contable.lg();
                
            }

            if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
            
            

        } 
    
    }
    
    
    
    
    LOG << "Finished Unit Test" << endl;
    
    return 0;
}
