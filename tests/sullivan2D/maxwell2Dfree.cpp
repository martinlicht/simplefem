

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
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.inclsullivan.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"
#include "../../vtk/vtkwriter.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: 2D Maxwell System" << nl;
    
    if(true){

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.getCoordinates().scale( 1.1 );
        
        M.check();
        
        M.automatic_dirichlet_flags();

        
        
        // const Float xfeq = 2.;
        // const Float yfeq = 2.;
        
        std::function<FloatVector(const FloatVector&)> experiment_sol = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 2 );
                // return FloatVector({ 1. });
                return FloatVector({ 
                    bumpfunction(vec[0])*bumpfunction(vec[1]) //std::sin( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
                    ,
                    bumpfunction(vec[0])*bumpfunction(vec[1]) //std::sin( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] )
                });
            };
        

        std::function<FloatVector(const FloatVector&)> experiment_ndiv = 
            [=](const FloatVector& vec) -> FloatVector {
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
            [=](const FloatVector& vec) -> FloatVector {
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
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 2 );
                
                

                
                
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
                    });
            };

            
            
            
        

        ConvergenceTable contable("Mass error and solver residual");
        
        contable << "sigma_error" << "u_error" << "du_error" << "residual" << "time" << nl;
        
        

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
            
            LOG << "Level: " << min_l << " <= " << l << " <= " << max_l << nl;
            LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            for( int r = min_r; r <= max_r; r++ )
            {
                
                LOG << "Level: "             << min_l << " <= " << l << " <= " << max_l << nl;
                LOG << "Polynomial degree: " << min_r << " <= " << r << " <= " << max_r << nl;
                        
                LOG << "... assemble matrices" << nl;
        
                
                auto scalar_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 ) );
                auto vector_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   ) );
                auto volume_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 ) );

                auto scalar_diffmatrix   = MatrixCSR( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 ) );
                auto scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                auto vector_diffmatrix   = MatrixCSR( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r ) );
                auto vector_diffmatrix_t = vector_diffmatrix.getTranspose();

                auto scalar_incmatrix   = MatrixCSR( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 ) );
                auto scalar_incmatrix_t = scalar_incmatrix.getTranspose();

                auto vector_incmatrix   = MatrixCSR( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   ) );
                auto vector_incmatrix_t = vector_incmatrix.getTranspose();

                auto volume_incmatrix   = MatrixCSR( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r-1 ) );
                auto volume_incmatrix_t = volume_incmatrix.getTranspose();
                
                auto mass = Conjugation( vector_massmatrix, vector_incmatrix );

                auto mat_A  = Conjugation( scalar_massmatrix, scalar_incmatrix );
                
                auto mat_Bt = scalar_incmatrix_t & scalar_diffmatrix_t & vector_massmatrix & vector_incmatrix; // upper right
                
                auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower left
                
                auto mat_C = Conjugation( volume_massmatrix, vector_diffmatrix & vector_incmatrix );
                
                LOG << "share zero A = " << mat_A.getnumberofzeroentries() << "/" <<  mat_A.getnumberofentries() << nl;
                LOG << "share zero B = " << mat_B.getnumberofzeroentries() << "/" <<  mat_B.getnumberofentries() << nl;
                LOG << "share zero C = " << mat_C.getnumberofzeroentries() << "/" <<  mat_C.getnumberofentries() << nl;
                
                auto& A  = mat_A;
                auto& Bt = mat_Bt;
                auto& B  = mat_B;
                auto& C  = mat_C;

                // auto negA  = A;  negA.scale(-1);
                // auto negB  = B;  negB.scale(-1);
                // auto negBt = Bt; negBt.scale(-1);
                
                auto SystemMatrix = C + B * inv(A,desired_precision) * Bt;
                
                {

                    const auto& function_ndiv = experiment_ndiv;
                    const auto& function_sol  = experiment_sol;
                    const auto& function_curl = experiment_curl;
                    const auto& function_rhs  = experiment_rhs;
                    
                    LOG << "... interpolate explicit solution and rhs" << nl;
                    
                    FloatVector interpol_ndiv = Interpolation( M, M.getinnerdimension(), 0, r+1, function_ndiv  );
                    FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 1, r,   function_sol  );
                    FloatVector interpol_curl = Interpolation( M, M.getinnerdimension(), 2, r-1, function_curl );
                    
                    FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 1, r,   function_rhs  );
                    
                    FloatVector rhs = vector_incmatrix_t * ( vector_massmatrix * interpol_rhs );

                        
                    FloatVector sol( vector_incmatrix.getdimin(), 0. );
                    
                    
                    timestamp start = timestampnow();
                    
                    {
                        
                        LOG << "... iterative solver" << nl;

                        auto PA = Conjugation( scalar_massmatrix, scalar_incmatrix )
                                  + 
                                  Conjugation( vector_massmatrix, scalar_diffmatrix & scalar_incmatrix );

                        auto PC = Conjugation( vector_massmatrix, vector_incmatrix )
                                  + 
                                  Conjugation( volume_massmatrix, vector_diffmatrix & vector_incmatrix );
                        
                        LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
                        LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;

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
                        
                        if(false){ // TODO(martinlicht): fix 
                    
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
                                desired_precision,
                                1,
                                desired_precision,
                                0
                            );
                        
                        }

                    }
                    
                    timestamp end = timestampnow();
    
                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                    
                    assert( sol.is_finite() );

                    auto ndiv = inv(A,desired_precision) * Bt * sol;
                    
                    auto curl = vector_diffmatrix * vector_incmatrix * sol;
                    
                    LOG << "... compute error and residual" << nl;

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

                    Float errornorm_ndiv = std::sqrt( errornorm_ndiv_sq );
                    Float errornorm_sol  = std::sqrt( errornorm_sol_sq  );
                    Float errornorm_curl = std::sqrt( errornorm_curl_sq );
                    
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

                    {
                        std::fstream fs( get_available_filename(get_basename(__FILE__)), std::fstream::out );
                        VTKWriter vtk( M, fs, get_basename(__FILE__) );

                        {
                            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r+1 );
                            const auto printable_sol = interpol_matrix * scalar_incmatrix * ndiv; 
                            vtk.write_cell_scalar_data( printable_sol, "computed_negative_divergence" , 1.0 );
                        }
                    
                        {
                            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r );
                            const auto printable_sol = interpol_matrix * vector_incmatrix * sol; 
                            vtk.write_cell_vector_data_barycentricgradients( printable_sol, "computed_solution" , 1.0 );
                        }
                    
                        {
                            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 2, 0, r-1 );
                            const auto printable_curl = interpol_matrix * vector_diffmatrix * vector_incmatrix * sol; 
                            Assert( printable_curl.getdimension() == (M.getinnerdimension()+1) * M.count_simplices(M.getinnerdimension()), 
                                            printable_curl.getdimension(), M.count_simplices(M.getinnerdimension()) );
                            vtk.write_cell_scalar_data_barycentricvolumes( printable_curl, "computed_curl" , 1.0 );
                        }
                        
                        vtk.write_cell_scalar_data( function_ndiv, "function_ndiv" , 1.0 );
                        vtk.write_cell_vector_data( function_sol,  "function_sol"  , 1.0 );
                        vtk.write_cell_scalar_data( function_curl, "function_curl" , 1.0 );

                        vtk.write_cell_vector_data( function_rhs,  "function_rhs" , 1.0 );
                    
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




