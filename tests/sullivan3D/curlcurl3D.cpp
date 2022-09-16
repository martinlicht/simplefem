

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
    
    LOG << "Unit Test for Solution of Darcy Problem" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = StandardCube3D();
    
    M.check();
    
    M.automatic_dirichlet_flags();
    M.check_dirichlet_flags();

    
    LOG << "Prepare scalar fields for testing..." << nl;
    

    const Float A = Constants::twopi;

            
    std::function<FloatVector(const FloatVector&)> experiment_sol = 
        [=](const FloatVector& vec) -> FloatVector{
            assert( vec.getdimension() == 3 );
            // return FloatVector({ 1. });
            return FloatVector({ 
                + blob( vec[0] ) * blob( vec[1] ) * blob_dev( vec[2] ) + blob( vec[0] ) * blob_dev( vec[1] ) * blob( vec[2] ) // u
                ,
                + blob( vec[0] ) * blob( vec[1] ) * blob_dev( vec[2] ) - blob_dev( vec[0] ) * blob( vec[1] ) * blob( vec[2] ) // v 
                ,
                - blob_dev( vec[0] ) * blob( vec[1] ) * blob( vec[2] ) - blob( vec[0] ) * blob_dev( vec[1] ) * blob( vec[2] ) // w 
                });
        };
    
    // u dx + v dy + w dz 
    // ->
    // u_y dyx + u_z dzx + v_x dxy + v_z dzy + w_x dxz + w_y dyz
    // =
    // - u_y dxy - u_z dxz + v_x dxy - v_z dyz + w_x dxz + w_y dyz
    // = 
    // ( v_x - u_y ) dxy + ( w_x - u_z ) dxz + ( w_y - v_z ) dyz
    
    std::function<FloatVector(const FloatVector&)> experiment_curl = 
        [=](const FloatVector& vec) -> FloatVector{
            assert( vec.getdimension() == 3 );
            // return FloatVector({ 1. });
            return FloatVector({ 
                + blob_dev( vec[0] ) * blob( vec[1] ) * blob_dev( vec[2] ) - blob_devdev( vec[0] ) * blob( vec[1] ) * blob( vec[2] ) //  v_x
                - blob( vec[0] ) * blob_dev( vec[1] ) * blob_dev( vec[2] ) - blob( vec[0] ) * blob_devdev( vec[1] ) * blob( vec[2] ) // -u_y
                ,
                - blob_devdev( vec[0] ) * blob( vec[1] ) * blob( vec[2] ) - blob_dev( vec[0] ) * blob_dev( vec[1] ) * blob( vec[2] ) //  w_x
                - blob( vec[0] ) * blob( vec[1] ) * blob_devdev( vec[2] ) - blob( vec[0] ) * blob_dev( vec[1] ) * blob_dev( vec[2] ) // -u_z
                ,
                - blob_dev( vec[0] ) * blob_dev( vec[1] ) * blob( vec[2] ) - blob( vec[0] ) * blob_devdev( vec[1] ) * blob( vec[2] ) //  w_y
                - blob( vec[0] ) * blob( vec[1] ) * blob_devdev( vec[2] ) + blob_dev( vec[0] ) * blob( vec[1] ) * blob_dev( vec[2] ) // -v_z
                });
        };
    
    // ( v_x - u_y ) dxy + ( w_x - u_z ) dxz + ( w_y - v_z ) dyz
    // -> (Hodge)
    // ( v_x - u_y ) dz  - ( w_x - u_z ) dy  + ( w_y - v_z ) dx
    // = 
    //  v_x dz - u_y dz  -  w_x dy + u_z dy  +  w_y dx - v_z dx
    // -> (Cartan)
    // v_xx dxz - u_xy dxz -  w_xx dxy + u_xz dxy  +  w_yy dyx - v_yz dyx
    // +
    // v_yx dyz - u_yy dyz -  w_xz dzy + u_zz dzy  +  w_zy dzx - v_zz dzx
    // = 
    // v_xx dxz - u_xy dxz -  w_xx dxy + u_xz dxy  -  w_yy dxy + v_yz dxy
    // +
    // v_yx dyz - u_yy dyz +  w_xz dyz - u_zz dyz  -  w_yy dxy + v_yz dxz
    // -> (Hodge)
    // - v_xx dy + u_xy dy -  w_xx dz + u_xz dz  -  w_yy dz + v_yz dz
    // +
    //   v_yx dx - u_yy dx +  w_xz dx - u_zz dx  +  w_yy dy - v_yz dy
    // = 
    // 
    //   v_yx dx - u_yy dx +  w_xz dx - u_zz dx
    // - v_xx dy + u_xy dy +  w_yy dy - v_yz dy
    // - w_xx dz + u_xz dz -  w_yy dz + v_yz dz
    // -> multiply by (-1)^(n(k-1)+1) = (-1)^(3+1)

    std::function<FloatVector(const FloatVector&)> experiment_rhs = 
        [=](const FloatVector& vec) -> FloatVector{
            assert( vec.getdimension() == 3 );
            return FloatVector({ 
                blob_dev(vec[0]) * blob(vec[1]) * blob(vec[2]) 
                -
                (
                + blob( vec[0] ) * blob( vec[1] ) * blob_dev( vec[2] ) + blob( vec[0] ) * blob_dev( vec[1] ) * blob( vec[2] ) // u
                + blob( vec[0] ) * blob( vec[1] ) * blob_dev( vec[2] ) + blob( vec[0] ) * blob_dev( vec[1] ) * blob( vec[2] ) // u
                + blob( vec[0] ) * blob( vec[1] ) * blob_dev( vec[2] ) + blob( vec[0] ) * blob_dev( vec[1] ) * blob( vec[2] ) // u
                ) 
                ,
                blob(vec[0]) * blob_dev(vec[1]) * blob(vec[2]) 
                + 
                // ( + blob_devdevdev( vec[0] ) * blob_dev( vec[1] ) + blob_dev( vec[0] ) * blob_devdev( vec[1] ) )
                (   24. * vec[0] * square( vec[1]*vec[1] - 1. )      + ( 4. * vec[0] * ( vec[0]*vec[0] - 1. ) ) * ( 12. * vec[1]*vec[1] - 4. ) )
                ,
                blob(vec[0]) * blob(vec[1]) * blob_dev(vec[2])
                + 
                // ( + blob_devdevdev( vec[0] ) * blob_dev( vec[1] ) + blob_dev( vec[0] ) * blob_devdev( vec[1] ) )
                (   24. * vec[0] * square( vec[1]*vec[1] - 1. )      + ( 4. * vec[0] * ( vec[0]*vec[0] - 1. ) ) * ( 12. * vec[1]*vec[1] - 4. ) )
                });
        };
    
    std::function<FloatVector(const FloatVector&)> experiment_aux = 
        [=](const FloatVector& vec) -> FloatVector{
            assert( vec.getdimension() == 3 );
            // return FloatVector({ 1. });
            return FloatVector( { 
                blob(vec[0])*blob(vec[1])*blob(vec[2])
                });
        };
    
    
    
    
    

    

    LOG << "Solving Poisson Problem with Neumann boundary conditions" << nl;

    const int min_l = 0; 
    const int max_l = 3;
    
    const int min_r = 1;
    const int max_r = 1;

    const int r_plus_vector = 1;
    const int r_plus_scalar = 1;
    
    
    ConvergenceTable contable("Mass error and numerical residuals");
    
    contable << "u_error" << "sigma_error" << "u_res" << "sigma_res" << "time";
    

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
                    
            LOG << "integration with: " << r_plus_vector << ", " << r_plus_scalar << nl;
                    
            LOG << "... assemble matrices" << nl;
    
            SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 + r_plus_scalar );
            SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   + r_plus_vector);
            SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

            SparseMatrix scalar_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

            SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

            SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

            SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
            SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

            SparseMatrix scalar_elevmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r+1, r_plus_scalar );
            SparseMatrix scalar_elevmatrix_t = scalar_elevmatrix.getTranspose();

            SparseMatrix vector_elevmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r  , r_plus_vector );
            SparseMatrix vector_elevmatrix_t = vector_elevmatrix.getTranspose();

            auto mat_A  = vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix;
            mat_A.sortandcompressentries();
                
            auto mat_Bt = vector_incmatrix_t & vector_elevmatrix_t & vector_massmatrix & vector_elevmatrix & scalar_diffmatrix & scalar_incmatrix; // upper right
            mat_Bt.sortandcompressentries();
            
            auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower bottom
            mat_B.sortandcompressentries();
            
            auto A  = MatrixCSR( mat_A  );
            auto Bt = MatrixCSR( mat_Bt );
            auto B  = MatrixCSR( mat_B  );
            
            auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
            
            // TODO: develop preconditioners 
            // auto PA = IdentityMatrix( A.getdimin() );
            // auto PC = IdentityMatrix( C.getdimin() );
            auto PA = MatrixCSR( vector_incmatrix_t & vector_elevmatrix_t & vector_massmatrix & vector_elevmatrix & vector_incmatrix )
                              + MatrixCSR( vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix );
            auto PC = MatrixCSR( scalar_incmatrix_t & scalar_diffmatrix_t & vector_elevmatrix_t & vector_massmatrix & vector_elevmatrix & scalar_diffmatrix & scalar_incmatrix );
                
            LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
            LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;
                        
                        
            const auto& function_sol = experiment_sol;
            const auto& function_rhs = experiment_rhs;
            const auto& function_aux = experiment_aux;
            
            LOG << "...interpolate explicit solution and rhs" << nl;
            
            FloatVector interpol_sol = Interpolation( M, M.getinnerdimension(), 1, r   + r_plus_vector, function_sol );
            FloatVector interpol_rhs = Interpolation( M, M.getinnerdimension(), 1, r   + r_plus_vector, function_rhs );
            FloatVector interpol_aux = Interpolation( M, M.getinnerdimension(), 0, r+1 + r_plus_scalar, function_aux );
            
            FloatVector rhs_sol = vector_incmatrix_t * vector_elevmatrix_t * vector_massmatrix * interpol_rhs;
            FloatVector rhs_aux = scalar_incmatrix_t * scalar_diffmatrix_t * vector_elevmatrix_t * vector_massmatrix * interpol_sol ;// FloatVector( B.getdimout(), 0. );
            
            FloatVector sol( A.getdimout(), 0. );
            FloatVector aux( B.getdimout(), 0. );

            // compute the solution ....



            timestamp start = gettimestamp();

            {
                const auto PAinv = inv(PA,desired_precision,-1);
                const auto PCinv = inv(PC,desired_precision,-1);
                BlockHerzogSoodhalterMethod( 
                    sol, 
                    aux, 
                    rhs_sol, 
                    rhs_aux, 
                    A, Bt, B, C, 
                    desired_precision,
                    1,
                    PAinv, PCinv
                );
            }

            timestamp end = gettimestamp();

            // ... computed the solution

            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
            
            LOG << "...compute error and residual:" << nl;

            auto errornorm_aux_sol = interpol_sol - vector_elevmatrix * vector_incmatrix * sol;
            auto errornorm_aux_aux = interpol_aux - scalar_elevmatrix * scalar_incmatrix * aux;

            Float errornorm_sol = sqrt( errornorm_aux_sol * ( vector_massmatrix * errornorm_aux_sol ) );
            Float errornorm_aux = sqrt( errornorm_aux_aux * ( scalar_massmatrix * errornorm_aux_aux ) );
            
            Float residual_sol  = ( rhs_sol - A * sol - Bt * aux ).norm();
            Float residual_aux  = ( rhs_aux - B * sol -  C * aux ).norm();

            LOG << "error:        " << errornorm_sol << nl;
            LOG << "aux error:    " << errornorm_aux << nl;
            LOG << "residual:     " << residual_sol << nl;
            LOG << "aux residual: " << residual_aux << nl;

            LOG << "aux rhs: " << rhs_aux.norm() << nl;

            contable << errornorm_sol;
            contable << errornorm_aux;
            contable << residual_sol;
            contable << residual_aux;
            contable << Float( end - start );
            contable << nl;

            contable.lg();
            
        }

        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
        
        

    } 

    
    LOG << "Finished Unit Test" << nl;
    
    return 0;
}
