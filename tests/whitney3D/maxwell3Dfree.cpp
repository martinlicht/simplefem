

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
// #include "../../solver/cgm.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/pcrm.hpp"
// #include "../../solver/minres.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: 3D Maxwell System" << nl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial3D M = StandardCube3D();
            
            M.getcoordinates().scale( 1.1 );
            
            M.check();
            
            M.automatic_dirichlet_flags();

            
            
            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [=](const FloatVector& vec) -> FloatVector{
                    Assert( vec.getdimension() == 3 );
                    // return FloatVector({ 1. });
                    return FloatVector({ 
                        bumpfunction(vec[0])*bumpfunction(vec[1])*bumpfunction(vec[2])
                        ,
                        bumpfunction(vec[0])*bumpfunction(vec[1])*bumpfunction(vec[2])
                        ,
                        bumpfunction(vec[0])*bumpfunction(vec[1])*bumpfunction(vec[2])
                    });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_ndiv = 
                [=](const FloatVector& vec) -> FloatVector{
                    Assert( vec.getdimension() == 3 );
                    // return FloatVector({ 1. });
                    return FloatVector( { 
                        - bumpfunction_dev(vec[0]) * bumpfunction(vec[1]) * bumpfunction(vec[2])
                        - bumpfunction(vec[0]) * bumpfunction_dev(vec[1]) * bumpfunction(vec[2])
                        - bumpfunction(vec[0]) * bumpfunction(vec[1]) * bumpfunction_dev(vec[2])
                    });
                };
            

            // + a_y dyx + a_z dzx + b_x dxy + b_z dzy + c_x dxz + c_y dyz
            // - a_y dxy - a_z dxz + b_x dxy - b_z dyz + c_x dxz + c_y dyz
            // - a_y dxy + b_x dxy - a_z dxz + c_x dxz - b_z dyz + c_y dyz
            // 

            std::function<FloatVector(const FloatVector&)> experiment_curl = 
                [=](const FloatVector& vec) -> FloatVector{
                    Assert( vec.getdimension() == 3 );
                    // return FloatVector({ 1. });
                    return FloatVector( { // - partial_y + partial_x
                        - bumpfunction(vec[0]) * bumpfunction_dev(vec[1]) * bumpfunction(vec[2]) // xy
                        + bumpfunction_dev(vec[0]) * bumpfunction(vec[1]) * bumpfunction(vec[2])
                        ,
                        - bumpfunction(vec[0]) * bumpfunction(vec[1]) * bumpfunction_dev(vec[2]) // xz
                        + bumpfunction_dev(vec[0]) * bumpfunction(vec[1]) * bumpfunction(vec[2])
                        ,
                        - bumpfunction(vec[0]) * bumpfunction(vec[1]) * bumpfunction_dev(vec[2]) // yz
                        + bumpfunction(vec[0]) * bumpfunction_dev(vec[1]) * bumpfunction(vec[2])
                    });
                };            

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    Assert( vec.getdimension() == 3 );
                    
                    return FloatVector({
                        -
                        bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1]) *        bumpfunction(vec[2])
                        -
                        bumpfunction(vec[0])        * bumpfunction_devdev(vec[1]) *        bumpfunction(vec[2])
                        -
                        bumpfunction(vec[0])        * bumpfunction(vec[1]) *        bumpfunction_devdev(vec[2])
                        ,
                        -
                        bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1]) *        bumpfunction(vec[2])
                        -
                        bumpfunction(vec[0])        * bumpfunction_devdev(vec[1]) *        bumpfunction(vec[2])
                        -
                        bumpfunction(vec[0])        * bumpfunction(vec[1]) *        bumpfunction_devdev(vec[2])
                        ,
                        -
                        bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1]) *        bumpfunction(vec[2])
                        -
                        bumpfunction(vec[0])        * bumpfunction_devdev(vec[1]) *        bumpfunction(vec[2])
                        -
                        bumpfunction(vec[0])        * bumpfunction(vec[1]) *        bumpfunction_devdev(vec[2])
                     });
                };

                
                
                
            

            ConvergenceTable contable("Mass error and solver residual");
            
            contable << "sigma_error" << "u_error" << "du_error" << "residual" << "residual" << "time";
            
            

            const int min_l = 0; 
            const int max_l = 4;
            
            const int min_r = 2; 
            const int max_r = 2;
            
            const int r_plus_scalar = 2;
            const int r_plus_vector = 2; 
            const int r_plus_pseudo = 2; 
            

            
            assert( 0 <= min_l and min_l <= max_l );
            assert( 0 <= min_r and min_r <= max_r );
            assert( 0 <= r_plus_scalar and 0 <= r_plus_vector and 0 <= r_plus_pseudo );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ )
            {
                
                LOG << "Level: " << l << "/" << max_l << nl;
                LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                for( int r = min_r; r <= max_r; r++ )
                {
                    
                    LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                    
                    //  LOG << "integration with: " << r_plus_scalar << ", " << r_plus_vector << ", " << r_plus_pseudo << nl;
                    
                    LOG << "... assemble matrices" << nl;
            
                    
                     // TODO: correct the degrees, perhaps via degree elevation
                    
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r + r_plus_scalar );
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r + r_plus_vector );
                    SparseMatrix pseudo_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r + r_plus_pseudo );

                    SparseMatrix scalar_augdiffmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r-1, 1 )
                                                          & FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );
                    SparseMatrix scalar_augdiffmatrix_t = scalar_augdiffmatrix.getTranspose();

                    SparseMatrix vector_augdiffmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r-1, 1 )
                                                          & FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix vector_augdiffmatrix_t = vector_augdiffmatrix.getTranspose();

                    SparseMatrix scalar_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r );
                    SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

                    SparseMatrix vector_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                    SparseMatrix pseudo_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 2, r );
                    SparseMatrix pseudo_incmatrix_t = pseudo_incmatrix.getTranspose();
                    
                    SparseMatrix scalar_elevmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_plus_scalar );
                    SparseMatrix scalar_elevmatrix_t = scalar_elevmatrix.getTranspose();

                    SparseMatrix vector_elevmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r, r_plus_vector );
                    SparseMatrix vector_elevmatrix_t = vector_elevmatrix.getTranspose();

                    SparseMatrix pseudo_elevmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r, r_plus_pseudo );
                    SparseMatrix pseudo_elevmatrix_t = pseudo_elevmatrix.getTranspose();

                

                    auto mat_A  = scalar_incmatrix_t & scalar_elevmatrix_t & scalar_massmatrix & scalar_elevmatrix & scalar_incmatrix;
                    mat_A.sortandcompressentries();
                    
                    auto mat_Bt = scalar_incmatrix_t & scalar_augdiffmatrix_t & vector_elevmatrix_t & vector_massmatrix & vector_elevmatrix & vector_incmatrix; // upper right
                    mat_Bt.sortandcompressentries();
                    
                    auto mat_B = mat_Bt.getTranspose(); //pseudo_incmatrix_t & pseudo_massmatrix & vector_elevmatrix & diffmatrix & vector_incmatrix; // lower bottom
                    mat_B.sortandcompressentries();
                    
                    auto mat_C  = vector_incmatrix_t & vector_augdiffmatrix_t & pseudo_elevmatrix_t & pseudo_massmatrix & pseudo_elevmatrix & vector_augdiffmatrix & vector_incmatrix;
                    mat_C.sortandcompressentries();
                    
                    LOG << "share zero A = " << mat_A.getnumberofzeroentries() << "/" << (Float) mat_A.getnumberofentries() << nl;
                    LOG << "share zero B = " << mat_B.getnumberofzeroentries() << "/" << (Float) mat_B.getnumberofentries() << nl;
                    LOG << "share zero C = " << mat_C.getnumberofzeroentries() << "/" << (Float) mat_C.getnumberofentries() << nl;
                    
                    auto A  = MatrixCSR( mat_A  );
                    auto Bt = MatrixCSR( mat_Bt );
                    auto B  = MatrixCSR( mat_B  );
                    auto C  = MatrixCSR( mat_C  );

                    auto negA  = -A;  
                    // auto negB  = -B;  
                    // auto negBt = -Bt; 
                    
                    {

                        const auto& function_ndiv = experiment_ndiv;
                        const auto& function_sol  = experiment_sol;
                        const auto& function_curl = experiment_curl;
                        const auto& function_rhs  = experiment_rhs;
                        
                        LOG << "...interpolate explicit solution and rhs" << nl;
                        
                        FloatVector interpol_ndiv = Interpolation( M, M.getinnerdimension(), 0, r + r_plus_scalar, function_ndiv );
                        FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 1, r + r_plus_vector, function_sol  );
                        FloatVector interpol_curl = Interpolation( M, M.getinnerdimension(), 2, r + r_plus_pseudo, function_curl );
                        
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 1, r + r_plus_vector, function_rhs  );
                        
                        FloatVector rhs = vector_incmatrix_t * vector_elevmatrix_t * ( vector_massmatrix * interpol_rhs );

                            
                        FloatVector sol( vector_incmatrix.getdimin(), 0. );
                        
                        
                        timestamp start = gettimestamp();
                        
                        // {

                            LOG << "...iterative solver" << nl;
                            
                            auto PA = MatrixCSR( scalar_incmatrix_t & scalar_elevmatrix_t & scalar_massmatrix & scalar_elevmatrix & scalar_incmatrix )
                                      + MatrixCSR( scalar_incmatrix_t & scalar_augdiffmatrix_t & vector_elevmatrix_t & vector_massmatrix & vector_elevmatrix & scalar_augdiffmatrix & scalar_incmatrix );
                            auto PC = MatrixCSR( vector_incmatrix_t & vector_elevmatrix_t & vector_massmatrix & vector_elevmatrix & vector_incmatrix )
                                      + MatrixCSR( vector_incmatrix_t & vector_augdiffmatrix_t & pseudo_elevmatrix_t & pseudo_massmatrix & pseudo_elevmatrix & vector_augdiffmatrix & vector_incmatrix );
                            
                            const auto PAinv = inv(PA,desired_precision,-1);
                            const auto PCinv = inv(PC,desired_precision,-1);

                            FloatVector  x_A( A.getdimin(),  0. ); 
                            FloatVector& x_C = sol;
                            
                            const FloatVector  b_A( A.getdimin(),  0. ); 
                            const FloatVector& b_C = rhs; 
                            
                            // auto Z  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix

                            BlockHerzogSoodhalterMethod( 
                                x_A, 
                                x_C, 
                                b_A, 
                                b_C, 
                                negA, Bt, B, C, 
                                desired_precision,
                                1,
                                PAinv, PCinv
                            );

                        // }
                        
                        timestamp end = gettimestamp();
        
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                        
                        assert( sol.isfinite() );

                        auto ndiv = x_A; // inv(A,1e-14) * Bt * sol;
                        
                        auto curl = vector_augdiffmatrix * vector_incmatrix * sol;
                        
                        LOG << "...compute error and residual:" << nl;

                        auto errornorm_aux_ndiv = interpol_ndiv - scalar_elevmatrix * scalar_incmatrix * ndiv;
                        auto errornorm_aux_sol  = interpol_sol  - vector_elevmatrix * vector_incmatrix *  sol;
                        auto errornorm_aux_curl = interpol_curl - pseudo_elevmatrix *                    curl;

                        Float errornorm_ndiv_sq = ( errornorm_aux_ndiv * ( scalar_massmatrix * errornorm_aux_ndiv ) );
                        Float errornorm_sol_sq  = ( errornorm_aux_sol  * ( vector_massmatrix * errornorm_aux_sol  ) );
                        Float errornorm_curl_sq = ( errornorm_aux_curl * ( pseudo_massmatrix * errornorm_aux_curl ) );
                        
                        LOG << errornorm_ndiv_sq << space << errornorm_sol_sq << space << errornorm_curl_sq << nl;

                        assert( errornorm_ndiv_sq >= 0 );
                        assert( errornorm_sol_sq  >= 0 );
                        assert( errornorm_curl_sq >= 0 );

                        Float errornorm_ndiv = sqrt( errornorm_ndiv_sq );
                        Float errornorm_sol  = sqrt( errornorm_sol_sq  );
                        Float errornorm_curl = sqrt( errornorm_curl_sq );
                        
                        LOG << "ndiv error sq: " << errornorm_ndiv_sq << nl;
                        LOG << "sol  error sq: " << errornorm_sol_sq << nl;
                        LOG << "curl error sq: " << errornorm_curl_sq << nl;
                        
                        LOG << "ndiv error: " << errornorm_ndiv << nl;
                        LOG << "sol  error: " << errornorm_sol << nl;
                        LOG << "curl error: " << errornorm_curl << nl;
                        
                        Float residualnorm_sol  = ( rhs - B * ndiv -  C * sol ).norm();
                        Float residualnorm_ndiv = (     - A * ndiv + Bt * sol ).norm();

                        LOG << "residual sol:  " << residualnorm_sol  << nl;
                        LOG << "residual ndiv: " << residualnorm_ndiv << nl;

                        contable << errornorm_ndiv;
                        contable << errornorm_sol;
                        contable << errornorm_curl;
                        
                        contable << residualnorm_sol;
                        contable << residualnorm_ndiv;
                        
                        contable << Float( end - start );
                        contable << nl;

                        contable.lg();

                    }
                    
                }

                if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }

            } 
            
            contable.lg();
        
        }
        
        
        
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}




