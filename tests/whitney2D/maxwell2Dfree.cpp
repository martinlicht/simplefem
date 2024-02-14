

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/math.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"
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

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: 2D Maxwell System" << nl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.getcoordinates().scale( 1.1 );
            
            M.check();
            
            M.automatic_dirichlet_flags();

            
            
            const Float xfeq = 2.;
            const Float yfeq = 2.;
            
            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [=](const FloatVector& vec) -> FloatVector{
                    Assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector({ 
                        bumpfunction(vec[0])*bumpfunction(vec[1]) //std::sin( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
                        ,
                        bumpfunction(vec[0])*bumpfunction(vec[1]) //std::sin( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
                    });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_ndiv = 
                [=](const FloatVector& vec) -> FloatVector{
                    Assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector( { 
                        - bumpfunction_dev(vec[0]) * bumpfunction(vec[1]) - bumpfunction(vec[0]) * bumpfunction_dev(vec[1])
//                             -xfeq * Constants::pi * std::cos( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
//                             +
//                             -yfeq * Constants::pi * std::sin( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
                    });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_curl = 
                [=](const FloatVector& vec) -> FloatVector{
                    Assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector( { // - partial_y + partial_x
                        - bumpfunction(vec[0]) * bumpfunction_dev(vec[1]) + bumpfunction_dev(vec[0]) * bumpfunction(vec[1])
//                             -
//                             yfeq * Constants::pi * 
//                             std::sin( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
//                             +
//                             xfeq * Constants::pi * 
//                             std::cos( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
                    });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    Assert( vec.getdimension() == 2 );
                    
                    
//                     const Float stepsize = desired_closeness;
//                     FloatVector ret(2);
//                     auto point = vec;
//                     
//                     FloatVector mid    = experiment_sol( point                              );
//                     FloatVector left   = experiment_sol( point - stepsize * unitvector(2,0) );
//                     FloatVector right  = experiment_sol( point + stepsize * unitvector(2,0) );
//                     FloatVector up     = experiment_sol( point - stepsize * unitvector(2,1) );
//                     FloatVector down   = experiment_sol( point + stepsize * unitvector(2,1) );
//                     
//                     ret = - ( up + down + left + right - 4 * mid );
//                     
//                     ret /= ( stepsize * stepsize );
//                     
//                     Assert( ret.getdimension() == 2 );
//                     
//                     return ret;

                    
                    
                    return FloatVector({
                        -
                        bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1])
                        -
                        bumpfunction(vec[0])        * bumpfunction_devdev(vec[1])
                        ,
                        -
                        bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1])
                        -
                        bumpfunction(vec[0])        * bumpfunction_devdev(vec[1])
                        
//                         ( xfeq*xfeq + yfeq*yfeq ) * Constants::pisquare 
//                         * 
//                         std::sin( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
//                         ,
//                         ( xfeq*xfeq + yfeq*yfeq ) * Constants::pisquare 
//                         * 
//                         std::sin( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
                     });
                };

                
                
                
            

            ConvergenceTable contable("Mass error and solver residual");
            
            contable << "sigma_error" << "u_error" << "du_error" << "residual" << "time";
            
            

            const int min_l = 0; 
            
            const int max_l = 4;
            
            const int min_r = 1; 
            
            const int max_r = 1;
            

            
            assert( 0 <= min_l and min_l <= max_l );
            assert( 0 <= min_r and min_r <= max_r );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ )
            {
                
                LOG << "Level: " << l << "/" << max_l << nl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                for( int r = min_r; r <= max_r; r++ )
                {
                    
                    LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                    
                    LOG << "... assemble matrices" << nl;
            
                    
                     // TODO: correct the degrees, perhaps via degree elevation
                    
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );

                    SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );
                    SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                    SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

                    SparseMatrix scalar_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r );
                    SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

                    SparseMatrix vector_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                    SparseMatrix vector_elevationmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r-1, 1 );
                    SparseMatrix vector_elevationmatrix_t = vector_elevationmatrix.getTranspose();

                    SparseMatrix volume_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 2, r );
                    SparseMatrix volume_incmatrix_t = volume_incmatrix.getTranspose();
                    
                    SparseMatrix volume_elevationmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r-1, 1 );
                    SparseMatrix volume_elevationmatrix_t = volume_elevationmatrix.getTranspose();

                

                    auto mass = vector_incmatrix_t * vector_massmatrix * vector_incmatrix;

                    auto mat_A  = scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix;
                    mat_A.sortandcompressentries();
                    
                    auto mat_Bt = scalar_incmatrix_t & scalar_diffmatrix_t & vector_elevationmatrix_t & vector_massmatrix & vector_incmatrix; // upper right
                    mat_Bt.sortandcompressentries();
                    
                    auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & vector_elevationmatrix & diffmatrix & vector_incmatrix; // lower bottom
                    mat_B.sortandcompressentries();
                    
                    auto mat_C  = vector_incmatrix_t & vector_diffmatrix_t & volume_elevationmatrix_t & volume_massmatrix & volume_elevationmatrix & vector_diffmatrix & vector_incmatrix;
                    mat_C.sortandcompressentries();
                    
                    LOG << "share zero A = " << mat_A.getnumberofzeroentries() << "/" <<  mat_A.getnumberofentries() << nl;
                    LOG << "share zero B = " << mat_B.getnumberofzeroentries() << "/" <<  mat_B.getnumberofentries() << nl;
                    LOG << "share zero C = " << mat_C.getnumberofzeroentries() << "/" <<  mat_C.getnumberofentries() << nl;
                    
                    auto A  = MatrixCSR( mat_A  );
                    auto Bt = MatrixCSR( mat_Bt );
                    auto B  = MatrixCSR( mat_B  );
                    auto C  = MatrixCSR( mat_C  );

                    // auto negA  = A;  negA.scale(-1);
                    // auto negB  = B;  negB.scale(-1);
                    // auto negBt = Bt; negBt.scale(-1);
                    
                    auto SystemMatrix = C + B * inv(A,desired_precision) * Bt;
                    
                    {

                        const auto& function_ndiv = experiment_ndiv;
                        const auto& function_sol  = experiment_sol;
                        const auto& function_curl = experiment_curl;
                        const auto& function_rhs  = experiment_rhs;
                        
                        LOG << "...interpolate explicit solution and rhs" << nl;
                        
                        FloatVector interpol_ndiv = Interpolation( M, M.getinnerdimension(), 0, r, function_ndiv );
                        FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 1, r, function_sol  );
                        FloatVector interpol_curl = Interpolation( M, M.getinnerdimension(), 2, r, function_curl );
                        
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 1, r, function_rhs  );
                        
                        FloatVector rhs = vector_incmatrix_t * ( vector_massmatrix * interpol_rhs );

                            
                        FloatVector sol( vector_incmatrix.getdimin(), 0. );
                        
                        
                        timestamp start = timestampnow();
                        
                        {

                            LOG << "...iterative solver" << nl;
                            
                            auto PA = MatrixCSR( scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix )
                                      + MatrixCSR( scalar_incmatrix_t & scalar_diffmatrix_t & vector_elevationmatrix_t & vector_massmatrix & vector_elevationmatrix & scalar_diffmatrix & scalar_incmatrix );
                            auto PC = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )
                                      + MatrixCSR( vector_incmatrix_t & vector_diffmatrix_t & vector_elevationmatrix_t & volume_massmatrix & vector_elevationmatrix & vector_diffmatrix & vector_incmatrix );
                            
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
                                -A, Bt, B, C, 
                                desired_precision,
                                1,
                                PAinv, PCinv
                            );

                        }
                        
                        timestamp end = timestampnow();
        
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                        
                        assert( sol.isfinite() );

                        auto ndiv = inv(A,desired_precision) * Bt * sol;
                        
                        auto curl = volume_elevationmatrix * vector_diffmatrix * vector_incmatrix * sol;
                        
                        LOG << "...compute error and residual:" << nl;

                        auto errornorm_aux_ndiv = interpol_ndiv - scalar_incmatrix * ndiv;
                        auto errornorm_aux_sol  = interpol_sol  - vector_incmatrix *  sol;
                        auto errornorm_aux_curl = interpol_curl -                    curl;

                        Float errornorm_ndiv_sq = ( errornorm_aux_ndiv * ( scalar_massmatrix * errornorm_aux_ndiv ) );
                        Float errornorm_sol_sq  = ( errornorm_aux_sol  * ( vector_massmatrix * errornorm_aux_sol  ) );
                        Float errornorm_curl_sq = ( errornorm_aux_curl * ( volume_massmatrix * errornorm_aux_curl ) );
                        Float residualnorm      = ( rhs - B * inv(A,desired_precision) * Bt * sol - C * sol ).norm();

                        LOG << errornorm_ndiv_sq << space << errornorm_sol_sq << space << errornorm_curl_sq << nl;

                        assert( errornorm_ndiv_sq >= 0 );
                        assert( errornorm_sol_sq  >= 0 );
                        assert( errornorm_curl_sq >= 0 );

                        Float errornorm_ndiv = sqrt( errornorm_ndiv_sq );
                        Float errornorm_sol  = sqrt( errornorm_sol_sq  );
                        Float errornorm_curl = sqrt( errornorm_curl_sq );
                        
                        LOG << "div  error sq: " << errornorm_ndiv_sq << nl;
                        LOG << "sol  error sq: " << errornorm_sol_sq << nl;
                        LOG << "curl error sq: " << errornorm_curl_sq << nl;
                        
                        LOG << "div  error: " << errornorm_ndiv << nl;
                        LOG << "sol  error: " << errornorm_sol << nl;
                        LOG << "curl error: " << errornorm_curl << nl;
                        
                        LOG << "residual:   " << residualnorm  << nl;

                        contable << errornorm_ndiv;
                        contable << errornorm_sol;
                        contable << errornorm_curl;
                        contable << residualnorm;
                        contable << Float( end - start );
                        contable << nl;

                        contable.lg();
                        


                    }
                    
                }

                if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }

                contable << nl;
                
                contable.lg();
        
            } 
            
            contable.lg();
        
        }
        
        
        
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}




